/*
create_maf.nf module

This module uses the vcf2maf.pl script to convert the vep-annotated
VCF files into MAF files.

vcf2maf.pl version: 1.6.17.
*/

process CREATE_MAF {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'medTime'
    publishDir "${params.output_dir}/mutation_calls/deepsomatic/raw_maf", mode: 'copy', pattern: "*deepsomatic*maf*"
    publishDir "${params.output_dir}/mutation_calls/MuSE/raw_maf", mode: 'copy', pattern: "*MuSE*maf*"

    input:
    tuple val(sample_id), path(vcf)
    path ref_fasta
    path ref_fasta_index
    path ref_dict

    output:
    tuple val(sample_id), path("*.maf")

    script:
    """
    # save the base name to change parameters based on variant caller
    BASE_NAME=\$(basename "${vcf}")

    if [[ "\$BASE_NAME" == *"deepsomatic"* ]]; then
        OUTPUT_NAME="${sample_id}.deepsomatic.vep.maf"
    elif [[ "\$BASE_NAME" == *"MuSE"* ]]; then
        OUTPUT_NAME="${sample_id}.MuSE.vep.maf"
    fi

    # DeepSomatic's REHEADER output has a single column (named TUMOR, see
    # reheader.nf) and no NORMAL column at all — use the tumor-only vcf2maf
    # invocation, not the matched-pair one, or --vcf-normal-id "NORMAL"
    # would point at a column that doesn't exist.
    if [[ "\$BASE_NAME" == *"deepsomatic"* ]]; then
        perl /app/vcf2maf.pl \
            --inhibit-vep \
            --input-vcf ${vcf} \
            --output-maf "\$OUTPUT_NAME" \
            --tumor-id ${sample_id} \
            --vcf-tumor-id "TUMOR" \
            --ref-fasta ${ref_fasta} \
            --ncbi-build hg38 \
            --maf-center NathansonLab

    # MuSE: matched tumor/normal pair, native NORMAL/TUMOR columns
    else
        perl /app/vcf2maf.pl \
            --inhibit-vep \
            --input-vcf ${vcf} \
            --output-maf "\$OUTPUT_NAME" \
            --tumor-id ${sample_id} \
            --normal-id "NORMAL" \
            --vcf-tumor-id "TUMOR" \
            --vcf-normal-id "NORMAL" \
            --ref-fasta ${ref_fasta} \
            --ncbi-build hg38 \
            --maf-center NathansonLab
    fi
    """
}
