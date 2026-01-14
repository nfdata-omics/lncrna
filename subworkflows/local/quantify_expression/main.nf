//
// Quantify gene/lncRNA expression using alignment-based or pseudo-alignment methods
//

include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts/main'
include { HTSEQ_COUNT           } from '../../../modules/nf-core/htseq/count/main'
include { biotypeInGtf          } from '../../../subworkflows/local/utils_nfcore_lncrna_pipeline'

workflow QUANTIFY_EXPRESSION {

    take:
    bam           // tuple val(meta), path(bam)      - Aligned BAM files (for alignment-based)
    bam_index     // tuple val(meta), path(bai)      - BAM indices
    gtf           // path(gtf)                       - Final annotation GTF
    counts_method        // val                      - quantification method

    main:
    ch_versions = Channel.empty()
    ch_counts = Channel.empty()
    ch_summary = Channel.empty()

    // Set default to featurecounts
    def method = counts_method ?: 'featurecounts'
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    //
    // Method 1: featureCounts (alignment-based, default)
    //
    if (!params.skip_qc && biotype && method == 'featurecounts') {
        // ch_bam_with_index = bam.join(bam_index, remainder: true)
        //     .map { meta, bam, bai ->
        //         [ meta, bam, bai ?: [] ]
        //     }

        gtf
            .map { biotypeInGtf(it, biotype) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        bam
            .combine(gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_counts = SUBREAD_FEATURECOUNTS.out.counts
//        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
    }

    //
    // Method 2: HTSeq (alignment-based)
    //
    else if (method == 'htseq') {
        ch_bam_with_index = bam.join(bam_index)

        HTSEQ_COUNT (
            ch_bam_with_index,
            gtf
        )
        ch_counts = HTSEQ_COUNT.out.txt
//        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions)
    }
    else {
        ch_counts = Channel.empty()
        error "Invalid quantification method: ${counts_method}. Must be 'featurecounts' or 'htseq'."
    }

    emit:
    counts   = ch_counts               // tuple val(meta), path(counts)
    summary  = ch_summary              // tuple val(meta), path(summary)
    versions = ch_versions.ifEmpty(null)
}
