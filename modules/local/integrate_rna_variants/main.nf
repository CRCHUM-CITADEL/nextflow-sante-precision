process INTEGRATE_RNA_VARIANTS {
    tag { sample_id }

    container params.container_r

    input:
        tuple val(sample_id), path(som_rna_vcf), path(som_dna_maf)

    output:
        tuple val(sample_id), path("${sample_id}.somatic_rna.maf")

    script:
    """
    gen_integrate_rna_variants.R \
        -d $som_dna_maf \
        -r $som_rna_vcf \
        -o ${sample_id}.somatic_rna.maf \
        --min_depth=3 --min_vaf=0.05
    """
}
