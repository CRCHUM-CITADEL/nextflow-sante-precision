process MERGE_EXPRESSION_FILES_TO_CBIOPORTAL {
    publishDir "${params.outdir}", mode: 'copy'

    container params.container_r

    input:
        path(tpm_file_list)

    output:
        path "data_expression.txt"


    // TODO: fix input_dir and pattern parameters
    script:
    """
    gen_merge_expression_files_to_cbioportal.R \
    --input_files ${tpm_file_list.join(',')} \
    --output_file data_expression.txt \
    --fill_missing 0 \
    --strict
    """
}
