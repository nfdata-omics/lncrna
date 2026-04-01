include { FILTER_KNOWN_LNCRNA       } from '../../../modules/local/filter_known_lncrna'
include { RENAME_LNCRNA             } from '../../../modules/local/rename_lncrna'
include { CLASSIFY_LNCRNA           } from '../../../modules/local/classify_lncrna'
include { MERGE_FINAL_ANNOTATION    } from '../../../modules/local/merge_final_annotation'
include { BLASTX                    } from '../../../modules/local/blastx'
include { FILTER_BLAST_RESULTS      } from '../../../modules/local/filter_blast_results'
include { GTF_TO_FASTA              } from '../../../modules/local/gtf_to_fasta'

workflow SUMMARY_AND_CLASSIFY_LNCRNA {

take:
    novel_lncrna_gtf                    // tuple val(meta), path(gtf) - Novel lncRNAs from COMBINE_PREDICTIONS (IDENTIFY_NOVEL_LNCRNA.out.final_lncrna_gtf)
    known_lncrna_gtf                    // path(known lncrnas gtf)- Known lncRNAs reference PREPARE_GENOME.out.known_lncrna_gtf
    protein_coding_gtf_novel            // path(gtf) - Protein-coding genes reference // GTF_FILTER_PROTEIN_CODING.out.protein_gtf
    // merge_reference_gtf              // path(gtf) - Full reference GTF used for final merge
    fasta                               // path(fasta) - Genome reference // PREPARE_GENOME.out.fasta
    blast_protein_db                    // tuple val(meta), path(db) - // PREPARE_GENOME.out.blast_protein_db
    lncrna_fasta                        //tuple val(meta), path(fasta) //PREPARE_GENOME.out.lncrna_fasta

main:
    ch_versions = Channel.empty()

//
// Filter out lncRNAs already known (keep only truly novel)
//
FILTER_KNOWN_LNCRNA (
    novel_lncrna_gtf,       // IDENTIFY_NOVEL_LNCRNA.out.final_lncrna_gtf
    known_lncrna_gtf        // PREPARE_GENOME.out.known_lncrna_gtf
    )
ch_truly_novel_gtf = FILTER_KNOWN_LNCRNA.out.gtf
    .map { meta, gtf -> [ [id: 'lncrnas_classification'], gtf ] }
//ch_versions = ch_versions.mix(FILTER_KNOWN_LNCRNA.out.versions)

//
// CONVERT: GTF to FASTA
//
GTF_TO_FASTA (
    ch_truly_novel_gtf,    // FILTER_KNOWN_LNCRNA.out.gtf
    fasta                  // PREPARE_GENOME.out.fasta
    )
ch_truly_novel_fasta = GTF_TO_FASTA.out.fasta
//ch_versions = ch_versions.mix(GTF_TO_FASTA.out.versions)

//
// MODULE: BLAST analysis
//
BLASTX (
    ch_truly_novel_fasta,       // GTF_TO_FASTA.out.fasta
    blast_protein_db.map { meta, db -> db }, // DB path
    params.blast_evalue ?: '1e-5'
    )
ch_blast_hits = BLASTX.out.blast_hits
//ch_versions = ch_versions.mix(BLASTX.out.versions)

//
// MODULE: Filter BLAST results
//
FILTER_BLAST_RESULTS (
    ch_truly_novel_fasta,       // GTF_TO_FASTA.out.fasta
    ch_blast_hits,              // BLASTX.out.blast_hits
    )
//ch_blast_filtered_fasta = FILTER_BLAST_RESULTS.out.filtered_fasta
//ch_versions = ch_versions.mix(FILTER_BLAST_RESULTS.out.versions)

//
// Rename lncRNAs based on closest protein-coding gene
//
// RENAME_LNCRNA (
//     ch_truly_novel_gtf,              // FILTER_KNOWN_LNCRNA.out.gtf
//     protein_coding_gtf_novel      // GTF_FILTER_PROTEIN_CODING.out.protein_gtf
//     )

// ch_novel_lncrna_gtf = RENAME_LNCRNA.out.gtf
//ch_versions = ch_versions.mix(RENAME_LNCRNA.out.versions)

//
// Classify lncRNAs by genomic location
//
CLASSIFY_LNCRNA (
    ch_truly_novel_gtf,                  // FILTER_KNOWN_LNCRNA.out.gtf
    protein_coding_gtf_novel      // GTF_FILTER_PROTEIN_CODING.out.protein_gtf
    )
ch_classification = CLASSIFY_LNCRNA.out.classification
ch_stats = CLASSIFY_LNCRNA.out.stats
//ch_versions = ch_versions.mix(CLASSIFY_LNCRNA.out.versions)

//
// Merge lncRNAs with protein-coding genes
//
MERGE_FINAL_ANNOTATION (
    ch_truly_novel_gtf,                  // FILTER_KNOWN_LNCRNA.out.gtf
    protein_coding_gtf_novel,     // GTF_FILTER_PROTEIN_CODING.out.protein_gtf    //merge_reference_gtf,
    fasta                            // PREPARE_GENOME.out.fasta
    )
ch_final_gtf = MERGE_FINAL_ANNOTATION.out.gtf
ch_final_fasta = MERGE_FINAL_ANNOTATION.out.fasta
ch_novel_lncrna_fasta = MERGE_FINAL_ANNOTATION.out.lncrna_fasta
ch_protein_fasta = MERGE_FINAL_ANNOTATION.out.protein_fasta
//ch_versions = ch_versions.mix(MERGE_FINAL_ANNOTATION.out.versions)

emit:
final_gtf              = ch_final_gtf              // tuple val(meta), path(gtf)  - Combined annotation
final_fasta            = ch_final_fasta            // tuple val(meta), path(fasta) - Combined sequences
novel_lncrna_gtf       = ch_truly_novel_gtf            // tuple val(meta), path(gtf)  - Renamed lncRNAs only //lncrna_gtf
novel_lncrna_fasta     = ch_novel_lncrna_fasta           // tuple val(meta), path(fasta) - lncRNA sequences //lncrna_fasta
lncrna_classification  = ch_classification         // tuple val(meta), path(txt)  - Classification by type
blast_hits             = ch_blast_hits             // found hits
//blast_filtered_fasta   = ch_blast_filtered_fasta   // lncRNA sequences with no protein similitudes
protein_fasta          = ch_protein_fasta          // tuple val(meta), path(fasta) - Protein sequences
lncrna_stats           = ch_stats
versions               = ch_versions.ifEmpty(null)
}
