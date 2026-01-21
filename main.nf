#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfdata-omics/lncrna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfdata-omics/lncrna
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LNCRNA                  } from './workflows/lncrna'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_lncrna_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_lncrna_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_lncrna_pipeline'
include { checkMaxContigSize      } from './subworkflows/local/utils_nfcore_lncrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta            = getGenomeAttribute('fasta')
params.additional_fasta = getGenomeAttribute('additional_fasta')
params.transcript_fasta = getGenomeAttribute('transcript_fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.gene_bed         = getGenomeAttribute('bed12')
params.bbsplit_index    = getGenomeAttribute('bbsplit')
params.sortmerna_index  = getGenomeAttribute('sortmerna')
params.star_index       = getGenomeAttribute('star')
params.hisat2_index     = getGenomeAttribute('hisat2')
params.salmon_index     = getGenomeAttribute('salmon')
params.kallisto_index   = getGenomeAttribute('kallisto')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFDATAOMICS_LNCRNA {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genome files
    //

    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.gencode,
        params.ncrna_fasta,
        params.validate_scaffolds,
        params.filter_lncrna_gtf,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.ribo_database_manifest,
        params.skip_gtf_filter,
        params.skip_bbsplit,
        params.skip_sortmerna,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.remove_ribo_rna
    )
    //ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { checkMaxContigSize(it) }
    }

    //
    // WORKFLOW: Run pipeline
    //
    ch_samplesheet = samplesheet

    LNCRNA (
        samplesheet,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.kallisto_index,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.rrna_fastas,                     // ch_ribo_db
        PREPARE_GENOME.out.sortmerna_index,
        PREPARE_GENOME.out.splicesites,
        PREPARE_GENOME.out.cds_fasta,                       // ch_cds_ref
        PREPARE_GENOME.out.lncrna_fasta,                    // ch_lncrna_ref
        PREPARE_GENOME.out.known_lncrna_gtf,                // ch_known_lncrna_gtf
        PREPARE_GENOME.out.lncrna_fasta,                    // ch_lncrna_fasta
        PREPARE_GENOME.out.blast_protein_db,                // blast_protein_database
        !params.sortmerna_index && params.remove_ribo_rna   // make_sortmerna_index
    )

    //ch_versions = ch_versions.mix(LNCRNA.out.versions)
    emit:
    multiqc_report = LNCRNA.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFDATAOMICS_LNCRNA (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFDATAOMICS_LNCRNA.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
