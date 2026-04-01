include { FILTER_TRANSCRIPTS_LENGTH } from '../../../modules/local/filter_transcripts_length'
include { FILTER_TRANSCRIPTS_EXONS  } from '../../../modules/local/filter_transcripts_exons'
include { EXTRACT_MRNA_SEQUENCES    } from '../../../modules/local/extract_mrna_sequences'
include { CPAT_BUILD_MODEL          } from '../../../modules/local/cpat_build_model'
include { CPAT                      } from '../../../modules/local/cpat'
include { PLEK_RUN                  } from '../../../modules/local/plek_run'
include { PLEK_PARSE                } from '../../../modules/local/plek_parse'
include { FEELNC_CODPOT_RUN         } from '../../../modules/local/feelnc_codpot_run'
include { FEELNC_CODPOT_PARSE       } from '../../../modules/local/feelnc_codpot_parse'
include { COMBINE_PREDICTIONS       } from '../../../modules/local/combine_predictions'


workflow IDENTIFY_NOVEL_LNCRNA {
    take:
        tmap                                  // gffcompare tracking file
        //ch_fasta                            // genome FASTA
        ch_lncrna_candidates_gtf              // filtered GTF (i, u, x class codes) //STRINGTIE_WORKFLOW.out.lncrna_candidates
        ch_lncrna_candidates_fa               // transcript sequences FASTA // STRINGTIE_WORKFLOW.out.lncrna_fasta
        ch_cds_fasta                          // CDS ref from prepare genome PREPARE_GENOME.out.cds_fasta
        ch_lncrna_fasta                       // lncRNA from prepare genome PREPARE_GENOME.out.lncrna_fasta
        ch_mrna_fasta                         // PREPARE_GENOME.out.mrna_fasta


    main:
        ch_versions = Channel.empty()

        //
        // FILTER: Remove short transcripts (< 200 nt)
        //
        FILTER_TRANSCRIPTS_LENGTH (
            ch_lncrna_candidates_gtf, //STRINGTIE_WORKFLOW.out.lncrna_candidates
            //ch_transcripts_fa
            ch_lncrna_candidates_fa //STRINGTIE_WORKFLOW.out.lncrna_fasta
        )
        ch_filtered_length_gtf = FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_gtf
        ch_filtered_length_fa  = FILTER_TRANSCRIPTS_LENGTH.out.filtered_length_fasta
//        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS_LENGTH.out.versions)


        //
        // FILTER: Remove single-exon transcripts
        //
        FILTER_TRANSCRIPTS_EXONS (
            ch_filtered_length_gtf,
            ch_filtered_length_fa
        )
        ch_filtered_exons_gtf = FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_gtf
        ch_filtered_exons_fa  = FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
//        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS_EXONS.out.versions)


        //
        // CPAT: Prepare or build models
        //
        ch_cpat_hexamer = Channel.empty()
        ch_cpat_logit = Channel.empty()

        if (!params.skip_cpat && params.cpat_hexamer && params.cpat_logit_model) {
            // Use provided models
            ch_cpat_hexamer = Channel.fromPath(params.cpat_hexamer)
            ch_cpat_logit = Channel.fromPath(params.cpat_logit_model)
        } else if (ch_cds_fasta && ch_lncrna_fasta) {
            // Build models
            CPAT_BUILD_MODEL (
                ch_cds_fasta, //coding sequences only
                ch_lncrna_fasta
            )
            ch_cpat_hexamer = CPAT_BUILD_MODEL.out.hexamer
            ch_cpat_logit = CPAT_BUILD_MODEL.out.logit_model
//            ch_versions = ch_versions.mix(CPAT_BUILD_MODEL.out.versions)
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
                ch_filtered_exons_fa, // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
                ch_cpat_hexamer,
                ch_cpat_logit
            )
            ch_cpat_results = CPAT.out.cpat_results
//            ch_versions = ch_versions.mix(CPAT.out.versions)
        }

        //
        // FEELnc: FlExible Extraction of LncRNAs
        //
        if (!params.skip_feelnc) {
            FEELNC_CODPOT_RUN (
                ch_filtered_exons_fa,   // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
                ch_mrna_fasta           // PREPARE_GENOME.out.mrna_fasta
            )

            FEELNC_CODPOT_PARSE (
                FEELNC_CODPOT_RUN.out.codpot_full,
                ch_filtered_exons_fa        // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
            )

            ch_feelnc_results = FEELNC_CODPOT_PARSE.out.feelnc_results
        }

        //
        // PLEK: Predictor of lncRNAs and mRNAs based on k-mer
        //
        if (!params.skip_plek) {
            PLEK_RUN (
                ch_filtered_exons_fa        // FILTER_TRANSCRIPTS_EXONS.out.filtered_exon_fasta
            )

            PLEK_PARSE (
                PLEK_RUN.out.plek_raw
            )
            ch_plek_results = PLEK_PARSE.out.plek_results
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
        ch_lncrna_pred_summary = COMBINE_PREDICTIONS.out.lncrna_pred_summary
        ch_lncrna_report    = COMBINE_PREDICTIONS.out.report
//        ch_versions = ch_versions.mix(COMBINE_PREDICTIONS.out.versions)


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
        lncrna_pred_summary = ch_lncrna_pred_summary     // lncrna prediction
        versions            = ch_versions
}
