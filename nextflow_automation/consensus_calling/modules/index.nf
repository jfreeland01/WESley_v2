/*
index.nf

This module compresses and indexes the sorted vcfs.

samtools version: 1.10.
*/

process INDEX {
    tag "${sample_id}"
    publishDir "${params.base_dir}/mutation_calls/consensus/raw-vcfs", mode: 'copy', pattern: "*vcf*"
    label 'lowCpu'
    label 'lowMem'
    label 'shortTime'

    input:
    tuple val(sample_id), path(deepsomatic_vcf), path(muse_vcf)

    output:
    tuple val(sample_id), path("*deepsomatic*.vcf.gz"), path("*deepsomatic*.vcf.gz.tbi"), path("*MuSE*.vcf.gz"), path("*MuSE*.vcf.gz.tbi")

    script:
    """
    # compress vcfs
    bgzip $deepsomatic_vcf
    bgzip $muse_vcf

    # index compressed vcfs
    tabix "${deepsomatic_vcf}.gz"
    tabix "${muse_vcf}.gz"
    """

}
