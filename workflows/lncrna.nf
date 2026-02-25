/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText; getStarPercentMapped } from '../subworkflows/local/utils_nfcore_lncrna_pipeline'
include { samplesheetToList                     } from 'plugin/nf-schema'

//
// SUBWORKFLOWS
//
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS  } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { FASTQ_ALIGN_STAR                      } from '../subworkflows/nf-core/fastq_align_star/main'
include { FASTQ_ALIGN_HISAT2                    } from '../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_STAR   } from '../subworkflows/nf-core/bam_dedup_umi'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_HISAT2 } from '../subworkflows/nf-core/bam_dedup_umi'
include { SUBREAD_FEATURECOUNTS                 } from '../modules/nf-core/subread/featurecounts/main'
include { HTSEQ_COUNT                           } from '../modules/nf-core/htseq/count/main'
include { STRINGTIE_WORKFLOW                    } from '../subworkflows/local/stringtie'
include { IDENTIFY_NOVEL_LNCRNA                 } from '../subworkflows/local/identify_novel_lncrna'
include { SUMMARY_AND_CLASSIFY_LNCRNA           } from '../subworkflows/local/summary_and_classify_lncrna'
include { EVALUATE_FINAL_LNCRNA                 } from '../subworkflows/local/evaluate_final_lncrna'
include { QUANTIFY_EXPRESSION                   } from '../subworkflows/local/quantify_expression'
include { GENERATE_COUNT_MATRIX                 } from '../modules/local/generate_count_matrix'
include { ANALYSIS_CIS_TRANS                    } from '../modules/local/analysis_cis_trans'
include { DIFFERENTIAL_EXPRESSION               } from '../modules/local/differential_expression'
include { QUANTIFY_PSEUDO_ALIGNMENT             } from '../subworkflows/nf-core/quantify_pseudo_alignment'
include { GENERATE_LNCRNA_REPORT                } from '../modules/local/generate_lncrna_report'
include { MERGE_FINAL_ANNOTATION                } from '../modules/local/merge_final_annotation'
include { CLASSIFY_LNCRNA                      } from '../modules/local/classify_lncrna'
include { GTF_FILTER_PROTEIN_CODING            } from '../modules/local/gtf_filter_protein_coding'
include { BAM_MARKDUPLICATES_PICARD            } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { FASTQC as FASTQC_TRIM                } from '../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LNCRNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_versions             // channel: [ path(versions.yml) ]
    ch_fasta                // channel: path(genome.fasta)
    ch_gtf                  // channel: path(genome.gtf)
    ch_fai                  // channel: path(genome.fai)
    ch_chrom_sizes          // channel: path(genome.sizes)
    ch_gene_bed             // channel: path(gene.bed)
    ch_transcript_fasta     // channel: path(transcript.fasta)
    ch_star_index           // channel: path(star/index/)
    ch_hisat2_index         // channel: path(hisat2/index/)
    ch_salmon_index         // channel: path(salmon/index/)
    ch_kallisto_index       // channel: [ meta, path(kallisto/index/) ]
    ch_bbsplit_index        // channel: path(bbsplit/index/)
    ch_ribo_db              // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index      // channel: path(sortmerna/index/)
    ch_splicesites          // channel: path(genome.splicesites.txt)
    ch_cds_ref              // channel: path(cds.gf)
    ch_lncrna_ref           // channel: path(lncrna.gtf)
    ch_known_lncrna_gtf     // channel: path(gtf_filtered_lncrna.gtf)
    lncrna_fasta
    blast_protein_database  // channel: path(lblast_protein_database)
    make_sortmerna_index    // boolean: Whether to create an index before running sortmerna

    main:

    //ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // ==========================================
    // SUBWORKFLOW: FASTQ preprocessing subworkflow
    // ==========================================

    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_samplesheet,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc,
        params.skip_trimming,
        params.skip_umi_extract,
        params.skip_linting,
        !salmon_index_available,
        !params.sortmerna_index,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        false,
        params.remove_ribo_rna,
        (params.ribo_removal_tool ?: 'sortmerna'),
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    //ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count
    if (!params.skip_fastqc && !params.skip_trimming && params.trimmer == 'trimgalore') {
        FASTQC_TRIM(
            FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect{ it[1] })
    }

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    // ==========================================
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    // ==========================================

    ch_genome_bam          = channel.empty()
    ch_genome_bam_index    = channel.empty()
    ch_star_log            = channel.empty()
    ch_unaligned_sequences = channel.empty()
    ch_transcriptome_bam   = channel.empty()

    if (!params.skip_alignment && params.aligner?.contains('star')) {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        FASTQ_ALIGN_STAR(
            FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads,
            ch_star_index.map { [ [:], it ] },
            ch_gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            ch_fasta.map { [ [:], it ] },
            ch_transcript_fasta.map { [ [:], it ] }
        )

        ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
        ch_transcriptome_bam       = FASTQ_ALIGN_STAR.out.bam_transcript
        ch_transcriptome_bai       = FASTQ_ALIGN_STAR.out.bai_transcript
        //ch_versions                = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        ch_multiqc_files = ch_multiqc_files
            .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_STAR.out.csi
        }
        //ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        //
        // Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {

            BAM_DEDUP_UMI_STAR(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta.map { [ [:], it ] },
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = BAM_DEDUP_UMI_STAR.out.bam
            ch_transcriptome_bam = BAM_DEDUP_UMI_STAR.out.transcriptome_bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_STAR.out.bai
            //ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_STAR.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_STAR.out.multiqc_files)

        } else {
            if (params.markduplicates) {
                BAM_MARKDUPLICATES_PICARD(
                    ch_genome_bam,
                    ch_fasta.map { [ [:], it ] },
                    ch_fai.map { [ [:], it ] }
                )

                ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
                ch_genome_bam_index = params.bam_csi_index ? BAM_MARKDUPLICATES_PICARD.out.csi : BAM_MARKDUPLICATES_PICARD.out.bai

                ch_multiqc_files = ch_multiqc_files
                    .mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]})
            } else {
                ch_multiqc_files = ch_multiqc_files
                    .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
                    .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
                    .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
            }
        }

    }

    // ==========================================
    // SUBWORKFLOW: Alignment with HISAT2
    // ==========================================
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { [ [:], it ] },
            ch_splicesites.map { [ [:], it ] },
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = FASTQ_ALIGN_HISAT2.out.bai
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        //ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        //
        // Remove duplicate reads from BAM file based on UMIs
        //

        if (params.with_umi) {

            BAM_DEDUP_UMI_HISAT2(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta.map { [ [:], it ] },
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = BAM_DEDUP_UMI_HISAT2.out.bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_HISAT2.out.bai
            //ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_HISAT2.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_HISAT2.out.multiqc_files)
        } else {
            if (params.markduplicates) {
                BAM_MARKDUPLICATES_PICARD(
                    ch_genome_bam,
                    ch_fasta.map { [ [:], it ] },
                    ch_fai.map { [ [:], it ] }
                )

                ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
                ch_genome_bam_index = params.bam_csi_index ? BAM_MARKDUPLICATES_PICARD.out.csi : BAM_MARKDUPLICATES_PICARD.out.bai

                ch_multiqc_files = ch_multiqc_files
                    .mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]})
                    .mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]})
            } else {
                ch_multiqc_files = ch_multiqc_files
                    .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{it[1]})
                    .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{it[1]})
                    .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{it[1]})
            }
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_log
            .map { meta, align_log -> [ meta ] + getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        // Save status for workflow summary
        ch_map_status = ch_percent_mapped
            .map {
                meta, mapped, pass ->
                    return [ meta.id, pass ]
            }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv'))
    }

    // ==========================================
    // Transcript assembly using Stringtie, merge GTFs, filter by classcode and create
    // FASTA of the filtered GTF
    // ==========================================
    def ch_stringtie_merged = Channel.empty()
    def ch_lncrna_candidates_gtf = Channel.empty()
    def ch_lncrna_candidates_fa = Channel.empty()
    def ch_gffcompare_annotated = Channel.empty()
    def tmap = Channel.empty()

    if (params.novel_lncrnas) {
        STRINGTIE_WORKFLOW (
            ch_genome_bam,
            ch_gtf,         // path(gtf) - PREPARE_GENOME.out.gtf
            ch_fasta,       // path(fasta) - PREPARE_GENOME.out.fasta
            ch_fai          // path(fai) - PREPARE_GENOME.out.fai
        )
        ch_stringtie_merged = STRINGTIE_WORKFLOW.out.stringtie_gtf_merged
        ch_lncrna_candidates_gtf = STRINGTIE_WORKFLOW.out.lncrna_candidates
        ch_lncrna_candidates_fa = STRINGTIE_WORKFLOW.out.lncrna_fasta
        ch_gffcompare_annotated = STRINGTIE_WORKFLOW.out.gffcompare_annotated
        tmap = STRINGTIE_WORKFLOW.out.tmap
        //ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)
    }

    // ==========================================
    // CODING POTENTIAL AND NOVEL LNCRNAS
    // ==========================================
    def ch_final_annotation
    def ch_lncrna_fasta
    def ch_lncrna_gtf
    def ch_protein_fasta
    def ch_lncrna_classification
    def ch_lncrna_stats
    def ch_cpat_hexamer
    def ch_cpat_logit
    def ch_cpat_plot
    def ch_classification_plot
    def ch_cpat_lncrna_validation
    def ch_cpat_coding_validation

    if (params.novel_lncrnas) {
        IDENTIFY_NOVEL_LNCRNA (
            tmap,
            ch_fasta,
            ch_lncrna_candidates_gtf,
            ch_lncrna_candidates_fa,
            ch_cds_ref,
            ch_lncrna_ref
        )

        validated_lncrnas_gtf   = IDENTIFY_NOVEL_LNCRNA.out.final_lncrna_gtf
        ch_cpat_hexamer         = IDENTIFY_NOVEL_LNCRNA.out.cpat_hexamer
        ch_cpat_logit           = IDENTIFY_NOVEL_LNCRNA.out.cpat_logit

        GTF_FILTER_PROTEIN_CODING (
            ch_gtf
        )
        def ch_protein_coding_gtf_novel = GTF_FILTER_PROTEIN_CODING.out.protein_gtf

        SUMMARY_AND_CLASSIFY_LNCRNA (
            validated_lncrnas_gtf,
            ch_known_lncrna_gtf,
            ch_protein_coding_gtf_novel,
            ch_gtf,
            ch_fasta,
            blast_protein_database,
            lncrna_fasta
        )

        ch_final_annotation      = SUMMARY_AND_CLASSIFY_LNCRNA.out.final_gtf
        ch_lncrna_fasta          = SUMMARY_AND_CLASSIFY_LNCRNA.out.lncrna_fasta
        ch_lncrna_gtf            = SUMMARY_AND_CLASSIFY_LNCRNA.out.lncrna_gtf
        ch_protein_fasta         = SUMMARY_AND_CLASSIFY_LNCRNA.out.protein_fasta
        ch_lncrna_classification = SUMMARY_AND_CLASSIFY_LNCRNA.out.lncrna_classification
        ch_lncrna_stats          = SUMMARY_AND_CLASSIFY_LNCRNA.out.lncrna_stats

        EVALUATE_FINAL_LNCRNA (
            ch_final_annotation,
            ch_lncrna_gtf,
            ch_lncrna_fasta,
            ch_protein_fasta,
            ch_lncrna_classification,
            ch_cpat_hexamer,
            ch_cpat_logit
        )

        ch_cpat_lncrna_validation = EVALUATE_FINAL_LNCRNA.out.cpat_lncrna_results
        ch_cpat_coding_validation = EVALUATE_FINAL_LNCRNA.out.cpat_coding_results
        ch_cpat_plot              = EVALUATE_FINAL_LNCRNA.out.cpat_plot
        ch_classification_plot    = EVALUATE_FINAL_LNCRNA.out.classification_plot
    }
    else {
        def ch_known_lncrna_tuple = ch_known_lncrna_gtf.map { gtf_file ->
            def meta = [ id: 'lncrnas_classification' ]
            [ meta, gtf_file ]
        }

        GTF_FILTER_PROTEIN_CODING (
            ch_gtf
        )
        def ch_protein_coding_gtf = GTF_FILTER_PROTEIN_CODING.out.protein_gtf

        CLASSIFY_LNCRNA (
            ch_known_lncrna_tuple,
            ch_protein_coding_gtf
        )
        ch_lncrna_classification = CLASSIFY_LNCRNA.out.classification
        ch_lncrna_stats          = CLASSIFY_LNCRNA.out.stats
        ch_final_annotation = ch_gtf.map { gtf_file ->
            def meta = [ id: 'reference' ]
            [ meta, gtf_file ]
        }
        ch_lncrna_fasta     = Channel.empty()
        ch_protein_fasta    = Channel.empty()
        ch_lncrna_gtf       = ch_known_lncrna_tuple

        ch_cpat_hexamer         = Channel.empty()
        ch_cpat_logit           = Channel.empty()
        ch_cpat_lncrna_validation = Channel.empty()
        ch_cpat_coding_validation = Channel.empty()
        ch_cpat_plot              = channel.value(file("$projectDir/assets/NO_FILE_cpat"))
        ch_classification_plot    = channel.value(file("$projectDir/assets/NO_FILE_cpat"))
    }

    // ==========================================
    // Quantification step (Featurecounts/Htseq)
    // ==========================================
    QUANTIFY_EXPRESSION (
        ch_genome_bam,                                      // BAMs originales
        ch_genome_bam_index,                                // Índices BAM
        ch_final_annotation.map { meta, gtf -> gtf },       // GTF con lncRNAs
        params.counts_method ?: 'featurecounts'     // method
    )
    ch_counts = QUANTIFY_EXPRESSION.out.counts
    //ch_versions = ch_versions.mix(QUANTIFY_EXPRESSION.out.versions)

    //
    // Count Matrix
    //
    GENERATE_COUNT_MATRIX (
        ch_counts.collect { meta, counts -> counts },
        ch_final_annotation.map { meta, gtf -> gtf },
        params.counts_method ?: 'featurecounts'
        )
    ch_alignment_matrix     = GENERATE_COUNT_MATRIX.out.matrix
    ch_alignment_lncrna_matrix = GENERATE_COUNT_MATRIX.out.lncrna_matrix
    ch_alignment_gene_info  = GENERATE_COUNT_MATRIX.out.gene_info

    //
    // MODULE: Cis/Trans Analysis
    //
    ANALYSIS_CIS_TRANS (
        ch_alignment_matrix.map { matrix -> [ ['id': 'all_samples'], matrix ] },
        ch_final_annotation.map { meta, gtf -> gtf }
    )
    ch_cis_results = ANALYSIS_CIS_TRANS.out.cis_results
    //ch_versions = ch_versions.mix(ANALYSIS_CIS_TRANS.out.versions)

    // ==========================================
    // Pseudoalignment and quantification with Salmon (OPCIONAL)
    // ==========================================
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {

        if (params.pseudo_aligner == 'salmon') {
            ch_pseudo_index = ch_salmon_index
        } else {
            ch_pseudo_index = ch_kallisto_index
        }

        QUANTIFY_PSEUDO_ALIGNMENT (
            channel.of([ [:], file(params.input, checkIfExists: true) ]),
            ch_strand_inferred_filtered_fastq,
            ch_pseudo_index.first(),
            ch_transcript_fasta.first(),
            ch_final_annotation.map { meta, gtf -> gtf }.first(),
            params.gtf_group_features ?: 'gene_id',
            params.gtf_extra_attributes ?: 'gene_name',
            params.pseudo_aligner ?: 'salmon',
            false,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen ?: 150,
            params.kallisto_quant_fraglen_sd ?: 20
        )

        ch_lncrna_gene_counts = QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled
        ch_lncrna_tpm = QUANTIFY_PSEUDO_ALIGNMENT.out.tpm_gene
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc.collect{it[1]})
        //ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.versions)
    }

    // ==========================================
    // Differential expression analysis (CONDICIONAL)
    // ==========================================
    if (params.design_file && !params.skip_differential_expression) {
        ch_design_file = Channel.fromPath(params.design_file, checkIfExists: true)
        // DE alignment-based
        DIFFERENTIAL_EXPRESSION (
            ch_alignment_lncrna_matrix.map { matrix -> [ ['id': 'alignment', 'contrast': params.contrast ?: '' ], matrix ] },
            //ch_alignment_gene_info,
            ch_design_file,
            params.contrast,
            params.lfc,
            params.fdr,
            params.min_count,
            params.min_total,
            params.correction_mode,
            params.batch ?: 'null',
            params.ruv_k,
            params.ruv_controls,
            params.top_genes,
            params.de_method ?: 'edgeR'
        )
        ch_de_results = DIFFERENTIAL_EXPRESSION.out.results_xlsx
        ch_de_normalized = DIFFERENTIAL_EXPRESSION.out.normalized
        // Plots
        ch_mds_png = DIFFERENTIAL_EXPRESSION.out.mds_png
        ch_pca_png = DIFFERENTIAL_EXPRESSION.out.pca_png
        ch_heatmap_global = DIFFERENTIAL_EXPRESSION.out.heatmap_png
        }

    //
    // MODULE: Generate final report
    //
    // ch_de_results_for_report = params.design_file && !params.skip_differential_expression ?
    //     ch_de_results.map { meta, results -> results } :
    //     channel.empty()
    // ch_de_plots_for_report = params.design_file && !params.skip_differential_expression ?
    //     ch_mds_png
    //         .mix(ch_pca_png)
    //         .mix(ch_heatmap_global)
    //         .collect() :
    //     channel.empty()
    // GENERATE_LNCRNA_REPORT (
    //     ch_de_results_for_report.ifEmpty(file("$projectDir/assets/NO_FILE_de_results")),
    //     ch_de_plots_for_report.ifEmpty(file("$projectDir/assets/NO_FILE_de_plots")),
    //     ch_lncrna_classification.map { meta, file -> file },
    //     ch_lncrna_stats.map { meta, file -> file },
    //     ch_cpat_plot.ifEmpty(file("$projectDir/assets/NO_FILE_cpat")),
    //     ch_alignment_matrix.ifEmpty(file("$projectDir/assets/NO_FILE_counts")),
    //     ch_final_annotation.map { meta, gtf -> gtf },
    //     params.novel_lncrnas ? ch_lncrna_gtf.map { meta, gtf -> gtf } : channel.value(file("$projectDir/assets/NO_FILE_rename"))
    // )
    // ch_final_report = GENERATE_LNCRNA_REPORT.out.report
    // ch_pipeline_summary = GENERATE_LNCRNA_REPORT.out.summary

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'lncrna_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
