/*
keep_nonsynonymous.nf module

This module filters synonymous mutations and keeps nonsynonymous mutation calls only
by referencing the 'nonsynonymous.txt' file.

Ubuntu version: 20.04.
*/

process KEEP_NONSYNONYMOUS {
    tag "${sample_id}"
    label 'lowCpu'
    label 'lowMem'
    label 'shortTime'

    input:
    tuple val(sample_id), path(maf_file)
    path nonsynonymous_list

    output:
    tuple val(sample_id), path("*vep.nonsynonymous.maf")

    script:
    """
    # save the base name to change parameters based on variant caller
    BASE_NAME=\$(basename "${maf_file}")

    if [[ "\$BASE_NAME" == *"deepsomatic"* ]]; then
        OUTPUT_NAME="${sample_id}.deepsomatic.vep.nonsynonymous.maf"
    elif [[ "\$BASE_NAME" == *"MuSE"* ]]; then
        OUTPUT_NAME="${sample_id}.MuSE.vep.nonsynonymous.maf"
    fi

    # filter synonymous mutations
    grep -v \
    -w \
    -F \
    -f \
    ${nonsynonymous_list} $maf_file > "\$OUTPUT_NAME"
    """
}
