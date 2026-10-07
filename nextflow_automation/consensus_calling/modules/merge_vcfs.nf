/*
merge_vcfs.nf

This module merges the two callers' union-passing VCFs (see intersect.nf)
into one file using GATK.

GATK version: 4.2.0.0.
*/

process MERGE_VCFS {
    tag "${sample_id}"
    label 'medCpu'
    label 'lowMem'
    label 'shortTime'

    input:
    tuple val(sample_id), path(deepsomatic_vcf), path(muse_vcf)

    output:
    tuple val(sample_id), path("${sample_id}*vcf.gz"), path("${sample_id}*.tbi")

    script:
    """
    gatk MergeVcfs \
    -I $deepsomatic_vcf \
    -I $muse_vcf \
    -O "${sample_id}.consensus.merged.vcf.gz"
    """
}
