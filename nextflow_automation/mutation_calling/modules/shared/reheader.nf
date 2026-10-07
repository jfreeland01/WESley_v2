/*
reheader.nf module

This module renames each caller's own VCF for consistent per-caller output
naming ahead of CREATE_MAF. With only DeepSomatic (single tumor-sample
column) and MuSE (native NORMAL/TUMOR two-column output) remaining as
callers, neither needs column reordering here — that TCGB-ID-based
reordering was only ever needed for Mutect2's raw output, which this
pipeline no longer runs. See CHANGES.md.

Note: this is the per-individual-caller reheader, used for each caller's
own MAF output. It is independent of consensus_calling's reheader.nf,
which harmonizes DeepSomatic's and MuSE's column conventions against each
other ahead of consensus merging — a different problem with a different
fix, see that module's docstring.

bcftools version: 1.10.
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
        cp "${vcf}" "${sample_id}.deepsomatic.reheader.vcf"
    elif [[ "\$BASE_NAME" == *"MuSE"* ]]; then
        cp "${vcf}" "${sample_id}.MuSE.reheader.vcf"
    fi
    """
}
