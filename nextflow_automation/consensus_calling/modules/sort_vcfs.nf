/*
sort_vcfs.nf module

This module sorts the records in VCF files according to the order
of the contigs in the header/sequence dictionary and then by coordinate.

GATK Version: 4.2.0.0
*/

process SORT_VCFS {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'medTime'

    input:
    tuple val(sample_id), path(deepsomatic_vcf), path(muse_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.deepsomatic.gatk-sort.vcf"), path("${sample_id}.muse.gatk-sort.vcf")

    script:
    """
    gatk SortVcf -I $deepsomatic_vcf -O "${sample_id}.deepsomatic.gatk-sort.vcf"
    gatk SortVcf -I $muse_vcf -O "${sample_id}.muse.gatk-sort.vcf"
    """
}
