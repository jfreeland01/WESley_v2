/*
INTERSECT.nf module

Builds the consensus call set as the UNION of DeepSomatic and MuSE calls
(a variant is kept if either caller found it) — not an intersection. This
is a deliberate change from the pipeline's previous 2-of-3 "majority vote"
rule. Benchmarking against SEQC2 (see CHANGES.md) found this simple union
beats every N-of-M voting rule tested, including the 3-caller version of
this same rule: Mutect2 and VarScan2's extra corroboration mostly added
false positives, not new true positives, once DeepSomatic was in the mix.

`bcftools isec -n +1` with 2 inputs means "keep sites present in 1 or more
of the inputs" — i.e. the union. (The old rule used `-n +2` across 3
inputs for "majority of 3.")

bcftools version: 1.10.2.
*/

process INTERSECT {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'medTime'

    input:
    tuple val(sample_id), path(deepsomatic_vcf), path(deepsomatic_index), path(muse_vcf), path(muse_index)

    output:
    tuple val(sample_id), path("**/0000.vcf"), path("**/0001.vcf")

    script:
    """
    bcftools isec -n +1 -p \
    $sample_id \
    $deepsomatic_vcf \
    $muse_vcf
    """
}
