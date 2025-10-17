process DRAGEN_FUSION_SV_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'


    tag { sample_id }

    container params.container_r

    input:
        tuple val(sample_id), path(dragen_fusion)

    output:
        path "${sample_id}.data_sv.txt"

    script:
    """
    gen_format_dragen_fusion.R \
        -i $dragen_fusion \
        -o ${sample_id}.data_sv.txt \
        -s $sample_id
    """
}
