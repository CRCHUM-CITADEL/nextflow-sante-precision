process GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    // use sample_id for logging
    tag { sample_id }

    container params.container_r

    input:
    tuple val(sample_id), path(somatic_cnv_vcf)
    path fold_changes_per_gene_cnv

    output:
    path "${sample_id}_data_cna_hg38.seg"

    script:
    """
    gen_cbioportal_converter.R \
      --vcf $somatic_cnv_vcf \
      --tsv $fold_changes_per_gene_cnv \
      --sample_id $sample_id \
      --output_dir .
    """
}
