process CONVERT_CPSR_TO_MAF {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    tag { sample_id }

    container params.container_r

    input:
        tuple val(sample_id), path(som_dna_rna_maf)
        path ger_dna_tsv_gz


    output:
        path "${sample_id}.somatic_rna_germline.maf"


    script:
    """
    zcat $ger_dna_tsv_gz > tmp.tsv
    head -1 tmp.tsv > tmp.germline.cpsr.tsv
    awk -F"\\t" '\$52=="Pathogenic" || \$52=="Likely_Pathogenic" || \$52=="VUS"' tmp.tsv >> tmp.germline.cpsr.tsv

    rm tmp.tsv # to reduce size of work dir

    gen_convert_cpsr_to_maf.R \
        tmp.germline.cpsr.tsv \
        $som_dna_rna_maf \
        ${sample_id}.somatic_rna_germline.maf
    """

}
