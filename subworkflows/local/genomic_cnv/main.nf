include { EXTRACT_GENE_CNV_FOLD_CHANGES } from '../../../modules/local/extract_gene_cnv_fold_changes'
include { GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL } from '../../../modules/local/gene_cnv_fold_changes_to_cbioportal'

workflow GENOMIC_CNV {
    take:
        cnv_vcf // tuple (sample_id, filepath)
        ensembl_annotations // gene annotation file (e.g. : "/lustre06/project/6079532/citadelomique/resources/genomes/grch38/annotations/ensembl/biomart_grch38_ensembl_113.tsv")

    main:

        fold_change_per_gene_cnv = EXTRACT_GENE_CNV_FOLD_CHANGES(
            cnv_vcf,
            ensembl_annotations
            )

        cbioportal_genomic_cnv_files = GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL(
            cnv_vcf,
            fold_change_per_gene_cnv
            )

        cbioportal_genomic_cnv_merged = cbioportal_genomic_cnv_files.collectFile( name : 'data_cna_hg38.seg', storeDir: "${params.outdir}", keepHeader : true, skip: 1)


    emit:
        cbioportal_genomic_cnv_merged

}
