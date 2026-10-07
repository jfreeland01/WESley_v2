/*
reheader.nf module

Harmonizes DeepSomatic's and MuSE's VCF column conventions so they can be
merged by consensus (bcftools isec / GATK MergeVcfs require identical
sample columns across inputs — the same requirement the old Mutect2 fix
addressed, see CHANGES.md fix #2).

DeepSomatic reports a single sample column named with the real tumor
sample ID (no NORMAL column — the normal sample is used internally to
classify PASS/GERMLINE/RefCall, not carried through). MuSE reports two
columns natively named NORMAL and TUMOR. To make them compatible, this
module reduces both to a single column literally named TUMOR:
DeepSomatic's one column is renamed from its real sample name to TUMOR;
MuSE's NORMAL column is dropped, keeping only its TUMOR column.

This means the final consensus VCF carries only the tumor's genotype at
each site, not the normal's — a deliberate simplification, not an
oversight: these are by definition sites classified as somatic (absent in
the matched normal), so the normal's genotype at them is not informative
to retain, and neither caller in this pipeline gives you real one anyway.

bcftools version: 1.10.2.
*/

process REHEADER {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'shortTime'

    input:
    tuple val(sample_id), path(deepsomatic_vcf), path(muse_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.deepsomatic.reheader.vcf"), path("${sample_id}.MuSE.reheader.vcf")

    script:
    """
    # DeepSomatic: single column, real sample name -> rename to TUMOR
    current_name=\$(grep "^#CHROM" "${deepsomatic_vcf}" | awk '{print \$10}')
    echo "\$current_name TUMOR" > deepsomatic-reheader.txt
    bcftools reheader "${deepsomatic_vcf}" \\
        -s deepsomatic-reheader.txt \\
        -o "${sample_id}.deepsomatic.reheader.vcf"

    # MuSE: drop NORMAL, keep only its native TUMOR column
    bcftools view -s TUMOR "${muse_vcf}" -o "${sample_id}.MuSE.reheader.vcf"
    """
}
