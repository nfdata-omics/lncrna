include { CPAT as RERUN_CPAT_LNCRNA        } from '../../../modules/local/cpat'
include { CPAT as RERUN_CPAT_CODING        } from '../../../modules/local/cpat'
include { GENERATE_FINAL_STATISTICS        } from '../../../modules/local/generate_final_statistics'

workflow EVALUATE_FINAL_LNCRNA {

take:
    final_gtf              // tuple val(meta), path(gtf)  - Combined annotation (lncRNAs + proteins)
    lncrna_gtf             // tuple val(meta), path(gtf)  - LncRNAs only GTF
    lncrna_fasta           // tuple val(meta), path(fa)   - LncRNA sequences
    protein_fasta          // tuple val(meta), path(fa)   - Protein-coding sequences
    lncrna_classification  // tuple val(meta), path(txt)  - Classification file
    cpat_hexamer           // path(tsv)                   - CPAT hexamer table
    cpat_logit_model       // path(RData)                 - CPAT logit model

main:
    ch_versions = Channel.empty()

//
// Re-run CPAT on final lncRNAs (validation)
//
RERUN_CPAT_LNCRNA (
    lncrna_fasta,           // lncrna_fasta
    cpat_hexamer,
    cpat_logit_model
)
ch_cpat_lncrna_results = RERUN_CPAT_LNCRNA.out.cpat_results
//ch_versions = ch_versions.mix(RERUN_CPAT_LNCRNA.out.versions)

//
// Re-run CPAT on protein-coding genes (control)
//
RERUN_CPAT_CODING (
    protein_fasta,
    cpat_hexamer,
    cpat_logit_model
)
ch_cpat_coding_results = RERUN_CPAT_CODING.out.cpat_results
//ch_versions = ch_versions.mix(RERUN_CPAT_CODING.out.versions)

//
// Generate final statistics and report
//
GENERATE_FINAL_STATISTICS (
    final_gtf,
    lncrna_gtf,
    ch_cpat_lncrna_results.map { meta, file -> file },
    ch_cpat_coding_results.map { meta, file -> file },
    lncrna_classification.map { meta, file -> file }
)
ch_final_stats_report = GENERATE_FINAL_STATISTICS.out.final_stats_report
ch_final_stats_summary = GENERATE_FINAL_STATISTICS.out.final_stats_summary
ch_cpat_plot = GENERATE_FINAL_STATISTICS.out.cpat_plot
ch_cpat_classification_plot = GENERATE_FINAL_STATISTICS.out.cpat_classification_plot
//ch_versions = ch_versions.mix(GENERATE_FINAL_STATISTICS.out.versions)

emit:
// CPAT results
cpat_lncrna_results    = ch_cpat_lncrna_results    // tuple val(meta), path(txt) - CPAT results for lncRNAs
cpat_coding_results    = ch_cpat_coding_results    // tuple val(meta), path(txt) - CPAT results for proteins

// Final reports and statistics
final_stats_report     = ch_final_stats_report          // path(html) - HTML report
final_stats_summary    = ch_final_stats_summary         // path(tsv)  - Summary statistics
cpat_plot              = ch_cpat_plot                   // path(png)  - CPAT comparison plot
classification_plot    = ch_cpat_classification_plot    // path(png)  - Classification barplot

// Versions
versions               = ch_versions.ifEmpty(null)
}
