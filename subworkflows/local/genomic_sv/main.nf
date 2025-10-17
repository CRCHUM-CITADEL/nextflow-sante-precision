// include modules
include { DRAGEN_FUSION_SV_TO_CBIOPORTAL } from '../../../modules/local/dragen_fusion_sv_to_cbioportal'


workflow GENOMIC_SV {
    take:
        sv_vcf

    main:
        cbioportal_genomic_sv_files = DRAGEN_FUSION_SV_TO_CBIOPORTAL(
            sv_vcf
        )

        cbioportal_genomic_sv_merged = cbioportal_genomic_sv_files.collectFile( name : 'data_sv.txt', storeDir: "${params.outdir}", keepHeader : true, skip: 1)

    emit:
        cbioportal_genomic_sv_merged
}
