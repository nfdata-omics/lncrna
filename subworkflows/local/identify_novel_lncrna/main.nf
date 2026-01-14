include { FILTER_TRANSCRIPTS_LENGTH } from '../../../modules/local/filter_transcripts_length'
include { FILTER_TRANSCRIPTS_EXONS  } from '../../../modules/local/filter_transcripts_exons'
include { CPAT_BUILD_MODEL          } from '../../../modules/local/cpat_build_model'
include { CPAT                      } from '../../../modules/local/cpat'
include { FEELNC_CODPOT             } from '../../../modules/local/feelnc_codpot'
include { PLEK                      } from '../../../modules/local/plek'
include { COMBINE_PREDICTIONS       } from '../../../modules/local/combine_predictions'


workflow IDENTIFY_NOVEL_LNCRNA {
    take:
        tmap                // gffcompare tracking file
        ch_fasta            // genome FASTA
        ch_gtf              // filtered GTF (i, u, x class codes)
        ch_transcripts_fa   // transcript sequences FASTA
        ch_cds_ref          // CDS ref from prepare genome
        ch_lncrna_ref       // lncRNA from prepare genome


    main:
        ch_versions = Channel.empty()

        //
        // FILTER: Remove short transcripts (< 200 nt)
        //
        FILTER_TRANSCRIPTS_LENGTH (
            ch_gtf,
            ch_transcripts_fa
        )
        ch_filtered_length_gtf = FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_gtf
        ch_filtered_length_fa  = FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_fasta
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS_LENGTH.out.versions)


        //
        // FILTER: Remove single-exon transcripts
        //
        FILTER_TRANSCRIPTS_EXONS (
            ch_filtered_length_gtf,
            ch_filtered_length_fa
        )
        ch_filtered_exons_gtf = FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_gtf
        ch_filtered_exons_fa  = FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS_EXONS.out.versions)


        //
        // CPAT: Prepare or build models
        //
        ch_cpat_hexamer = Channel.empty()
        ch_cpat_logit = Channel.empty()

        if (!params.skip_cpat && params.cpat_hexamer && params.cpat_logit_model) {
            // Use provided models
            ch_cpat_hexamer = Channel.fromPath(params.cpat_hexamer)
            ch_cpat_logit = Channel.fromPath(params.cpat_logit_model)
        } else if (ch_cds_ref && ch_lncrna_ref) {
            // Build models
            CPAT_BUILD_MODEL (
                ch_cds_ref,
                ch_lncrna_ref
            )
            ch_cpat_hexamer = CPAT_BUILD_MODEL.out.hexamer
            ch_cpat_logit = CPAT_BUILD_MODEL.out.logit_model
            ch_versions = ch_versions.mix(CPAT_BUILD_MODEL.out.versions)
        } else {
            error "ERROR: Must provide either:\n" +
                    "  - Pre-built CPAT models (--cpat_hexamer and --cpat_logit_model)\n" +
                    "  OR\n" +
                    "  - Training data (coding_fasta and noncoding_fasta)"
        }

        //
        // CPAT: Coding Potential Assessment Tool
        //
        if (!params.skip_cpat) {
            CPAT (
                ch_filtered_exons_fa,
                ch_cpat_hexamer,
                ch_cpat_logit
            )
            ch_cpat_results = CPAT.out.cpat_results
            ch_versions = ch_versions.mix(CPAT.out.versions)
        }

        //
        // FEELnc: FlExible Extraction of LncRNAs
        //
        if (!params.skip_feelnc) {
            FEELNC_CODPOT (
                ch_filtered_exons_fa,
                ch_fasta
            )
            ch_feelnc_results = FEELNC_CODPOT.out.feelnc_results
            ch_versions = ch_versions.mix(FEELNC_CODPOT.out.versions)
        }


        //
        // PLEK: Predictor of lncRNAs and mRNAs based on k-mer
        //
        if (!params.skip_plek) {
            PLEK (
                ch_filtered_exons_fa
            )
            ch_plek_results = PLEK.out.plek_results
            ch_versions = ch_versions.mix(PLEK.out.versions)
        }

        //
        // COMBINE: Merge predictions from CPAT, FEELnc and PLEK
        //
        COMBINE_PREDICTIONS (
            ch_cpat_results,
            ch_feelnc_results,
            ch_plek_results,
            ch_filtered_exons_gtf,
            ch_filtered_exons_fa,
            tmap
        )
        ch_final_lncrna_gtf = COMBINE_PREDICTIONS.out.lncrna_gtf
        ch_final_lncrna_fa  = COMBINE_PREDICTIONS.out.lncrna_fasta
        ch_lncrna_report    = COMBINE_PREDICTIONS.out.report
        ch_versions = ch_versions.mix(COMBINE_PREDICTIONS.out.versions)


    emit:
        filtered_length_gtf = ch_filtered_length_gtf     // After length filter
        filtered_exons_gtf  = ch_filtered_exons_gtf      // After exon filter
        cpat_hexamer        = ch_cpat_hexamer            // CPAT hexamer (provided or built)
        cpat_logit          = ch_cpat_logit              // CPAT logit (provided or built)
        cpat_results        = ch_cpat_results            // CPAT predictions
        feelnc_results      = ch_feelnc_results          // FEELnc predictions
        plek_results        = ch_plek_results            // PLEK predictions
        final_lncrna_gtf    = ch_final_lncrna_gtf        // Final lncRNA GTF
        final_lncrna_fasta  = ch_final_lncrna_fa         // Final lncRNA FASTA
        report              = ch_lncrna_report           // Summary report
        versions            = ch_versions
}
