nextflow.enable.dsl=2

// import modules
include { MUSE } from './modules/muse/muse.nf'
include { DEEPSOMATIC } from './modules/deepsomatic/deepsomatic.nf'
include { INDEX } from './modules/shared/index.nf'
include { SELECT_VARIANTS } from './modules/shared/select_variants.nf'
include { VEP; VEP_OMICS } from './modules/shared/vep.nf'
include { REHEADER } from './modules/shared/reheader.nf'
include { CREATE_MAF } from './modules/shared/create_maf.nf'
include { KEEP_NONSYNONYMOUS } from './modules/shared/keep_nonsynonymous.nf'
include { RENAME_HG38 } from './modules/shared/rename_hg38.nf'
include { ONCOKB; ONCOKB_OMICS } from './modules/shared/oncokb.nf'
include { KEEP_TERTP } from './modules/shared/keep_tertP.nf'


// define functions for user guidance
def help_message() {
    log.info """Usage:

    The typical command for running the pipeline is as follows:

    nextflow -C <CONFIG_PATH> run mutation_calling.nf --output_dir <PATH> --ref_dir <PATH> -params-file manifest.json --capture_bed <PATH> [OPTIONS]

    Required arguments:
    --output_dir                  Path to the output directory for results
    --ref_dir                     Path to the reference directory
    --samples                     Path to manifest JSON file produced by make_mc_manifest.py
    --capture_bed                 Path to a plain BED file of exome capture regions

    Optional arguments:
    --cpus                        Number of CPUs to use for processing (default: 30)
    --help                        Show this help message and exit

    Generating the manifest (local):
    python make_mc_manifest.py --platform local --bam_dir /path/to/bams --metadata metadata.xlsx --output manifest.json

    Generating the manifest (HealthOmics):
    python make_mc_manifest.py --platform omics --store_id <ID> --region <REGION> --metadata metadata.xlsx --output manifest.json

    Examples:

    # Basic usage with required parameters
    nextflow -C nextflow.config \\
        run mutation_calling.nf \\
        --output_dir /path/to/data \\
        --ref_dir /path/to/reference \\
        --samples manifest.json \\
        --capture_bed /path/to/capture_regions.bed
    """.stripIndent()
}

def log_workflow() {
    log.info """\
 __     __     ______     ______     __         ______     __  __
/\\ \\  _ \\ \\   /\\  ___\\   /\\  ___\\   /\\ \\       /\\  ___\\   /\\ \\_\\ \\
\\ \\ \\/ ".\\ \\  \\ \\  __\\   \\ \\___  \\  \\ \\ \\____  \\ \\  __\\   \\ \\____ \\
 \\ \\__/".~\\_\\  \\ \\_____\\  \\/\\_____\\  \\ \\_____\\  \\ \\_____\\  \\/\\_____\\
  \\/_/   \\/_/   \\/_____/   \\/_____/   \\/_____/   \\/_____/   \\/_____/
=========================================================================================
    Workflow ran:       : ${workflow.manifest.name}
    Command ran         : ${workflow.commandLine}
    Started on          : ${workflow.start}
    Config File used    : ${workflow.configFiles ?: 'None specified'}
    Container(s)        : ${workflow.containerEngine}:${workflow.container ?: 'None'}
    Nextflow Version    : ${workflow.manifest.nextflowVersion}
    """.stripIndent()
}


def parameter_validation() {
    // Parameter validation
    if (!params.output_dir) {
        error "ERROR: --output_dir parameter is required"
        exit 1
    }

    if (!params.capture_bed) {
        error "ERROR: --capture_bed parameter is required"
        exit 1
    }
}


// main workflow
workflow {
    // Show help message if requested
    if (params.help) {
        help_message()
        exit(0)
    }
    // validate parameters
    parameter_validation()
    if (!params.samples) {
        error "ERROR: --samples is required. Generate a manifest with make_mc_manifest.py and pass it via --samples manifest.json"
        exit 1
    }

    // stage reference files from params (derive from ref_dir when not set explicitly)
    ref_fasta           = file(params.ref_fasta             ?: "${params.ref_dir}/Homo_sapiens_assembly38.fasta")
    ref_fasta_index     = file(params.ref_fasta_index       ?: "${params.ref_dir}/Homo_sapiens_assembly38.fasta.fai")
    ref_dict            = file(params.ref_dict              ?: "${params.ref_dir}/Homo_sapiens_assembly38.dict")
    muse_dbsnp          = file(params.muse_dbsnp            ?: "${params.ref_dir}/common_all_20180418.vcf.gz")
    muse_dbsnp_index    = file(params.muse_dbsnp_index      ?: "${params.ref_dir}/common_all_20180418.vcf.gz.tbi")
    capture_bed         = file(params.capture_bed)
    nonsynonymous_list  = file(params.nonsynonymous_list    ?: "${params.ref_dir}/nonsynonymous.txt")
    // vep_cache: on HealthOmics it's an S3 URI to be staged (path); locally
    // it's an in-container path string (val) so Nextflow doesn't bind-mount
    // /opt/vep over the container's VEP install. See VEP / VEP_OMICS in vep.nf.
    is_omics            = System.getenv('AWS_WORKFLOW_RUN') as Boolean
    vep_cache           = is_omics ? file(params.vep_cache) : params.vep_cache
    sample_list         = file(params.samples)

    // logging workflow details
    log_workflow()

    // channel in samples from manifest JSON file
    Channel.fromPath(sample_list)
        .splitJson(path: 'samples')
        .map { s ->
            tuple(
                s.sample_id,
                s.tumor_id,
                file(s.tumor_bam),
                s.tumor_bai  != null ? file(s.tumor_bai)  : [],
                s.tumor_sbi  != null ? file(s.tumor_sbi)  : [],
                s.normal_id  != null ? s.normal_id         : 'NO_FILE',
                s.normal_bam != null ? file(s.normal_bam)  : [],
                s.normal_bai != null ? file(s.normal_bai)  : []
            )
        }
        .set { bams }

    // split channels based on if normal is available — both remaining
    // callers (DeepSomatic, MuSE) are tumor/normal-paired only
    bams
    .branch { row ->
        paired: row[6] != []
        tumor_only: row[6] == []
    }
    .set { samples }

    // run MuSE variant caller
    muse_vcfs = MUSE(samples.paired, ref_fasta, ref_fasta_index, ref_dict, muse_dbsnp, muse_dbsnp_index)

    // run DeepSomatic variant caller
    deepsomatic_vcfs = DEEPSOMATIC(samples.paired, ref_fasta, ref_fasta_index, ref_dict, capture_bed)

    // concatenate both callers' data channels
    filtered_vcfs = deepsomatic_vcfs.concat(muse_vcfs)

    // compress and index the vcfs
    compressed_vcfs = INDEX(filtered_vcfs)

    // select for passing variants via gatk SelectVariants
    selected_vcfs = SELECT_VARIANTS(compressed_vcfs)

    // annotate for biological effects via VEP — dispatch to the path-input
    // variant on HealthOmics (S3 staging) or the val-input variant locally
    // (uses cache baked into e10m/vep:115 at /opt/vep/.vep)
    vep_annotated_vcfs = is_omics
        ? VEP_OMICS(selected_vcfs, ref_fasta, ref_fasta_index, ref_dict, vep_cache)
        : VEP(selected_vcfs, ref_fasta, ref_fasta_index, ref_dict, vep_cache)

    // change the column names in the vcf for standardization
    reheadered_vcfs = REHEADER(vep_annotated_vcfs)

    // generate MAF files
    maf_files = CREATE_MAF(reheadered_vcfs, ref_fasta, ref_fasta_index, ref_dict)

    // grep and keep only tert promoter mutations
    KEEP_TERTP(maf_files)

    // filter out synonymous mutations
    nonsynonymous_mutations = KEEP_NONSYNONYMOUS(maf_files, nonsynonymous_list)

    // rename and reformat the files
    renamed_files = RENAME_HG38(nonsynonymous_mutations)

    // oncokb annotation for clinical relevance — dispatch to the AWS Secrets
    // Manager variant on HealthOmics, otherwise use the Nextflow secret variant
    is_omics
        ? ONCOKB_OMICS(renamed_files)
        : ONCOKB(renamed_files)
}
