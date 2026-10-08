// keep_tertP.nf
//
// This module greps the raw MAF files for tert promoter (tertP) mutations.
//
// Ubuntu version: 20.04

process KEEP_TERTP {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'shortTime'
    publishDir "${params.output_dir}/mutation_calls/deepsomatic/tertp", mode: 'copy', pattern: "*deepsomatic*maf*"
    publishDir "${params.output_dir}/mutation_calls/MuSE/tertp", mode: 'copy', pattern: "*MuSE*maf*"

    input:
    tuple val(sample_id), path(maf)

    output:
    tuple val(sample_id), path("*tertp*.maf")

    script:
    // parse base name in Groovy
    def basename = file(maf).baseName

    // define output names
    def output_file = null

    if (basename.contains("deepsomatic")) {
        output_file = "${sample_id}.deepsomatic.vep.tertp.maf"
    } else if (basename.contains("MuSE")) {
        output_file = "${sample_id}.MuSE.vep.tertp.maf"
    }

    """
    # maintain header
    grep "^#\\|^Hugo_Symbol" $maf > $output_file || true

    # keep only TERT promoter mutations
    grep "TERT" $maf | grep "upstream_gene_variant" >> $output_file || true
    """
}
