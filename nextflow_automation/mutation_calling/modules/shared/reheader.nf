/*
reheader.nf module

This module renames each caller's own VCF for consistent per-caller output
naming ahead of CREATE_MAF. DeepSomatic reports a single sample column
named with the real tumor sample ID (no NORMAL column — see
deepsomatic.nf's docstring); this module renames that one column to the
literal name TUMOR so downstream modules (CREATE_MAF in particular, which
reads `--vcf-tumor-id`) have one consistent column name to expect across
callers rather than branching on dynamic real sample names. MuSE already
reports a column literally named TUMOR, so it just gets copied through.

Note: this is the per-individual-caller reheader, used for each caller's
own MAF output. It is independent of consensus_calling's reheader.nf,
which does the analogous harmonization ahead of consensus merging instead.

bcftools version: 1.10.2.
*/

process REHEADER {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'shortTime'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}*reheader.vcf")

    script:
    """
    BASE_NAME=\$(basename "${vcf}")

    if [[ "\$BASE_NAME" == *"deepsomatic"* ]]; then
        current_name=\$(grep "^#CHROM" "${vcf}" | awk '{print \$10}')
        echo "\$current_name TUMOR" > deepsomatic-reheader.txt
        bcftools reheader "${vcf}" \\
            -s deepsomatic-reheader.txt \\
            -o "${sample_id}.deepsomatic.reheader.vcf"
    elif [[ "\$BASE_NAME" == *"MuSE"* ]]; then
        cp "${vcf}" "${sample_id}.MuSE.reheader.vcf"
    fi
    """
}
