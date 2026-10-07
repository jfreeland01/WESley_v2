/*
deepsomatic.nf module

Runs Google's DeepSomatic (deep-learning tumor/normal somatic variant caller).
Output is a plain (uncompressed) VCF named without the typical .vcf.gz
extension, which DeepSomatic's postprocess_variants step respects by writing
plain text — this keeps the module's output compatible with the shared
INDEX module (bgzip/tabix), exactly like MuSE's and VarScan2's raw output
did, rather than double-compressing.

Raw output includes RefCall/GERMLINE/NoCall records alongside PASS
(DeepSomatic reports on every examined site, not just variants). PASS-only
filtering happens later via the existing shared SELECT_VARIANTS step
(--exclude-filtered), the same place every other caller's output gets
filtered — see select_variants.nf.

DeepSomatic's tumor/normal mode reports a single sample column (the tumor's
real sample ID, via --sample_name_tumor) — unlike Mutect2/MuSE/VarScan2, it
does not emit a NORMAL genotype column (the normal sample is used
internally to classify PASS/GERMLINE/RefCall, not carried through as its
own column). See reheader.nf and CHANGES.md for how this is reconciled
with MuSE's native two-column NORMAL/TUMOR output before consensus merging.

DeepSomatic Version: 1.10.0
*/

process DEEPSOMATIC {
    tag "${sample_id}"
    label 'highCpu'
    label 'highMem'
    label 'extraLongTime'

    input:
    tuple val(sample_id), val(tumor_id), path(tumor_bam), path(tumor_bai), path(tumor_sbi), val(normal_id), path(normal_bam), path(normal_bai)
    path ref_fasta
    path ref_fasta_index
    path ref_dict
    path capture_bed

    output:
    tuple val(sample_id), val(tumor_id), val(normal_id), path("${sample_id}.deepsomatic.vcf")

    script:
    """
    run_deepsomatic \\
        --model_type=WES \\
        --ref=${ref_fasta} \\
        --reads_normal=${normal_bam} \\
        --reads_tumor=${tumor_bam} \\
        --output_vcf="${sample_id}.deepsomatic.vcf" \\
        --regions=${capture_bed} \\
        --num_shards=${task.cpus} \\
        --sample_name_normal="${normal_id}" \\
        --sample_name_tumor="${tumor_id}"
    """
}
