include { VCF2MAF } from '../../../modules/nf-core/vcf2maf'
include { INTEGRATE_RNA_VARIANTS } from '../../../modules/local/integrate_rna_variants'
include { PCGR } from '../../../modules/local/pcgr'
include { CONVERT_CPSR_TO_MAF } from '../../../modules/local/convert_cpsr_to_maf'
include { DOWNLOAD_VEP_TEST } from '../../../modules/local/download_vep_test'
include { DOWNLOAD_PCGR } from '../../../modules/local/download_pcgr'

workflow GENOMIC_MUTATIONS {
    take:
        ger_dna_vcf // tuple (sample_id, filepath)
        som_dna_vcf // tuple (sample_id, filepath)
        som_rna_vcf // tuple (sample_id, filepath)
        fasta
        vep_cache
        pcgr_data

    main:

        DOWNLOAD_VEP_TEST()
        DOWNLOAD_PCGR()

        // Combine both channels and take the first available
        ch_vep_data = vep_cache
            .mix(DOWNLOAD_VEP_TEST.out.cache_dir)
            .first()

        ch_pcgr_data = pcgr_data
            .mix(DOWNLOAD_PCGR.out.data_dir)
            .first()

        som_dna_vcf_input = som_dna_vcf.map { id, vcf ->
            def meta = [ id: id ]
            return tuple(meta, vcf)
        }

        VCF2MAF(
            som_dna_vcf_input,
            fasta,
            ch_vep_data,
            params.vep_params
        )

        som_dna_maf = VCF2MAF.out.maf.map { meta, vcf ->
            return tuple(meta.id, vcf)
        }

        // join on ID to create tuple(id, dna, rna)
        som_rna_dna_tuple = som_rna_vcf
            .join(som_dna_maf)

        som_dna_rna_maf = INTEGRATE_RNA_VARIANTS(
            som_rna_dna_tuple
        )

        ger_dna_tsv = PCGR (
            ger_dna_vcf,
            ch_vep_data,
            ch_pcgr_data
        )

        cbioportal_genomic_mutation_files = CONVERT_CPSR_TO_MAF (
            som_dna_rna_maf,
            ger_dna_tsv
        )

        cbioportal_genomic_mutations_merged = cbioportal_genomic_mutation_files.collectFile( name : 'data_mutations_dna_rna_germline.txt', storeDir: "${params.outdir}", keepHeader : true, skip: 1)


    emit:
        cbioportal_genomic_mutations_merged

}
