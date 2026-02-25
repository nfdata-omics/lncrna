//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_NCRNA_FASTA      } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
//include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx/main'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { GFFREAD as FASTA_EXTRACT_TRANSCRIPTS } from '../../../modules/nf-core/gffread'
include { GFFREAD as PROTEIN_EXTRACT        } from '../../../modules/nf-core/gffread'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { KALLISTO_INDEX                    } from '../../../modules/nf-core/kallisto/index'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'
include { GTF_FILTER_LNCRNA                    } from '../../../modules/local/gtf_filter_lncrna'
include { EXTRACT_CDS_SEQUENCES                } from '../../../modules/local/extract_cds_sequences'
include { EXTRACT_MRNA_SEQUENCES               } from '../../../modules/local/extract_mrna_sequences'
include { EXTRACT_LNCRNA_SEQUENCES             } from '../../../modules/local/extract_lncrna_sequences'

include { BLAST_MAKEBLASTDB                    } from '../../../modules/nf-core/blast/makeblastdb/main'

workflow PREPARE_GENOME {
    take:
    fasta                    //      file: /path/to/genome.fasta
    gtf                      //      file: /path/to/genome.gtf
    gff                      //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    gene_bed                 //      file: /path/to/gene.bed
    splicesites              //      file: /path/to/splicesites.txt
    bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    //sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    salmon_index             // directory: /path/to/salmon/index/
    kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    ncrna_fasta              //      file: /path/to/ncrnas.fasta [OPTIONAL]
    validate_scaffolds       //   boolean: Validate scaffolds exist in FASTA (default: true)
    filter_lncrna_gtf        //   boolean: Filter GTF to only lncRNA (default: false)
    featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'
    pseudo_aligner           //    string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    ribo_database_manifest  // directory: /path/to/ribo_database_manifest
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             //   boolean: Skip BBSplit for removal of non-reference genome reads
    skip_sortmerna           //   boolean: Skip sortmerna for removal of reads mapping to sequences in sortmerna_fasta_list
    skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    //   boolean: Skip all of the pseudoalignment-based processes within the pipeline
    remove_ribo_rna

    main:
    ch_versions         = Channel.empty()
    ch_known_lncrna_gtf = Channel.empty()

    //
    // Check if there is a FASTA
    //
    def fasta_provided = (fasta ? true : false)

    ch_fasta = Channel.of([])
    if (fasta_provided) {
    // Uncompress FASTA if needed
        if (fasta.endsWith('.gz')) {
            ch_fasta    = GUNZIP_FASTA ( [ [:], file(fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
//            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = Channel.value(file(fasta, checkIfExists: true))
        }
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], file(gtf, checkIfExists: true) ] ).gunzip.map { it[1] }
//                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf, checkIfExists: true))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], file(gff, checkIfExists: true) ] ).gunzip
//                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff, checkIfExists: true)).map { [ [:], it ] }
            }
            ch_gtf      = GFFREAD ( ch_gff, [] ).gtf.map { it[1] }
//            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        //
        // Filter GTF once (combines nf-core and lncRNA logic)
        //
        def should_filter_gtf =
            validate_scaffolds ||
            (
                ((!skip_alignment && aligner) ||
                (!skip_pseudo_alignment && pseudo_aligner) ||
                (!transcript_fasta)
                ) && (!skip_gtf_filter)
            )

        if (should_filter_gtf) {
            GTF_FILTER ( ch_fasta, ch_gtf )
            ch_gtf_validated = GTF_FILTER.out.genome_gtf
//            ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        } else {
            ch_gtf_validated = ch_gtf
        }

        // Filter to lncRNA only (optional)
        if (filter_lncrna_gtf) {
            GTF_FILTER_LNCRNA ( ch_gtf_validated )
            ch_gtf_to_use = GTF_FILTER_LNCRNA.out.lncrna_gtf
            ch_known_lncrna_gtf = GTF_FILTER_LNCRNA.out.lncrna_gtf
//            ch_versions = ch_versions.mix(GTF_FILTER_LNCRNA.out.versions)
        } else {
            ch_gtf_to_use = ch_gtf_validated
        }

        // Update main GTF reference
        ch_gtf = ch_gtf_validated
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : featurecounts_group_type
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], file(additional_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
//            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta, checkIfExists: true))
        }

        CUSTOM_CATADDITIONALFASTA (
            ch_fasta.combine(ch_gtf).map { fasta, gtf -> [ [:], fasta, gtf ] },
            ch_add_fasta.map { [ [:], it ] },
            biotype
        )
        ch_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta.map { it[1] }.first()
        ch_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf.map { it[1] }.first()
//        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], file(gene_bed, checkIfExists: true) ] ).gunzip.map { it[1] }
//            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed, checkIfExists: true))
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
//        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //----------------------------------------------------------------------
    // Transcript FASTA:
    //    - If provided, decompress (optionally preprocess if GENCODE)
    //    - If not provided but have genome+GTF, create from them
    //----------------------------------------------------------------------
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], file(transcript_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
//            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta, checkIfExists: true))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
//            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        // Extract transcripts from genome only if annotation is available
        if (gtf || gff) {
            ch_gtf_fasta_inputs = ch_gtf_to_use
                .combine(ch_fasta)
                .map { gtf, fasta -> [ [id:'transcripts'], gtf, fasta] }
                .multiMap { meta, gtf, fasta ->
                    gtf_input: [meta, gtf]
                    fasta_input: fasta
                }

            FASTA_EXTRACT_TRANSCRIPTS (
                ch_gtf_fasta_inputs.gtf_input,
                ch_gtf_fasta_inputs.fasta_input
            )
            ch_transcript_fasta = FASTA_EXTRACT_TRANSCRIPTS.out.gffread_fasta.map { meta, fasta ->
                if (fasta instanceof List) {
                    def expected = (meta?.id ? "${meta.id}.fasta" : null)
                    def chosen = expected ? fasta.find { it.getName() == expected } : fasta.find { it.getName().endsWith('.fasta') }
                    return chosen ?: fasta[0]
                }
                return fasta
            }
            // ch_versions = ch_versions.mix(FASTA_EXTRACT_TRANSCRIPTS.out.versions)
        } else {
            ch_transcript_fasta = Channel.empty()
        }
    }

    //
    // Concatenate with ncrna_fasta if provided
    //
    ch_ncrna_fasta = Channel.empty()
    if (ncrna_fasta) {
        if (ncrna_fasta.endsWith('.gz')) {
            ch_ncrna_fasta = GUNZIP_NCRNA_FASTA ( [ [:], file(ncrna_fasta, checkIfExists: true) ] ).gunzip.map { it[1] }
//            ch_versions    = ch_versions.mix(GUNZIP_NCRNA_FASTA.out.versions)
        } else {
            ch_ncrna_fasta = Channel.value(file(ncrna_fasta, checkIfExists: true))
        }

        ch_transcript_fasta = ch_transcript_fasta.combine(ch_ncrna_fasta)
            .map { transcripts, ncrna ->
                def outPath = "${transcripts.parent}/${transcripts.baseName}_with_ncrna.fasta"
                def combined_file = new File(outPath)
                combined_file.text = transcripts.text + "\n" + ncrna.text
                return file(outPath, checkIfExists: true)
            }
    }

    //
    // Create chromosome sizes file
    //
    // CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    // ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    // ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    // ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    SAMTOOLS_FAIDX (
    ch_fasta.map { [ [:], it ] },  // [meta, fasta]
    [[],[]],                       // [meta2, fai] - empty for new index
    true                           // get_sizes = true (generates .sizes file)
    )
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { it[1] }
//    ch_versions    = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // Get list of indices that need to be created
    //
    def prepare_tool_indices = []
    if (!skip_bbsplit)                              { prepare_tool_indices << 'bbsplit' }
    if (!skip_sortmerna)                            { prepare_tool_indices << 'sortmerna' }
    if (!skip_alignment) {
        if (aligner && aligner.contains('star')) {
            prepare_tool_indices << 'star'
        } else if (aligner == 'hisat2') {
            prepare_tool_indices << 'hisat2'
        }
    }
    def have_transcripts = (transcript_fasta || gtf || gff)
    if (!skip_pseudo_alignment && pseudo_aligner && have_transcripts)   { prepare_tool_indices << pseudo_aligner }

    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ [:], bbsplit_index ] ).untar.map { it[1] }
//                ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index))
            }
        } else if (fasta_provided && bbsplit_fasta_list) {
            Channel
                .from(file(bbsplit_fasta_list, checkIfExists: true))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fafile -> [ [ 'id', id ], [ 'fasta', file(fafile, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { it -> it[1] } // Get rid of keys and keep grouped values
                .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT( [ [:], [] ], [], ch_fasta, ch_bbsplit_fasta_list, true ).index
//            ch_versions      = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
        }
    }

    //
    // Uncompress sortmerna index or generate from scratch if required
    //
    ch_sortmerna_index = Channel.empty()
    ch_rrna_fastas = Channel.empty()

    if ('sortmerna' in prepare_tool_indices) {
        ribo_db = file(ribo_database_manifest)

        // SortMeRNA needs the rRNAs even if we're providing the index
        ch_rrna_fastas = Channel.from(ribo_db.readLines())
            .map { row -> file(row, checkIfExists: true) }

        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ( [ [:], sortmerna_index ] ).untar.map { it[1] }
//                ch_versions = ch_versions.mix(UNTAR_SORTMERNA_INDEX.out.versions)
            } else {
                ch_sortmerna_index = Channel.value([[:], file(sortmerna_index)])
            }
        } else {

            SORTMERNA_INDEX (
                Channel.of([ [],[] ]),
                ch_rrna_fastas.collect().map { [ 'rrna_refs', it ] },
                Channel.of([ [],[] ])
            )
            ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
//            ch_versions = ch_versions.mix(SORTMERNA_INDEX.out.versions)
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
//                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index))
            }
        } else {
            // Check if an AWS iGenome has been provided to use the appropriate version of STAR
            def is_aws_igenome = false
            if (fasta && gtf) {
                if ((file(fasta).getName() - '.gz' == 'genome.fa') && (file(gtf).getName() - '.gz' == 'genes.gtf')) {
                    is_aws_igenome = true
                }
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
//                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
//                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
//            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.value(file(splicesites))
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], hisat2_index ] ).untar.map { it[1] }
//                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(hisat2_index))
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
//            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
//            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
//            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    //
    // Uncompress Kallisto index or generate from scratch if required
    //
    ch_kallisto_index = Channel.empty()
    if (kallisto_index) {
        if (kallisto_index.endsWith('.tar.gz')) {
            ch_kallisto_index = UNTAR_KALLISTO_INDEX ( [ [:], kallisto_index ] ).untar
//            ch_versions     = ch_versions.mix(UNTAR_KALLISTO_INDEX.out.versions)
        } else {
            ch_kallisto_index = Channel.value([[:], file(kallisto_index)])
        }
    } else {
        if ('kallisto' in prepare_tool_indices) {
            ch_kallisto_index = KALLISTO_INDEX ( ch_transcript_fasta.map { [ [:], it] } ).index
//            ch_versions     = ch_versions.mix(KALLISTO_INDEX.out.versions)
        }
    }

    //
    // Extract reference sequences and build protein database
    //
    EXTRACT_CDS_SEQUENCES ( ch_fasta, ch_gtf )
    ch_cds_fasta = EXTRACT_CDS_SEQUENCES.out.fasta
//    ch_versions = ch_versions.mix(EXTRACT_CDS_SEQUENCES.out.versions)

    EXTRACT_MRNA_SEQUENCES ( ch_fasta, ch_gtf )
    ch_mrna_fasta = EXTRACT_MRNA_SEQUENCES.out.fasta
//    ch_versions = ch_versions.mix(EXTRACT_MRNA_SEQUENCES.out.versions)

    EXTRACT_LNCRNA_SEQUENCES ( ch_fasta, ch_gtf )
    ch_lncrna_fasta = EXTRACT_LNCRNA_SEQUENCES.out.fasta
//    ch_versions = ch_versions.mix(EXTRACT_LNCRNA_SEQUENCES.out.versions)

    ch_gtf_fasta_inputs_prot = ch_gtf
        .combine(ch_fasta)
        .map { gtf, fasta -> [ [id:'proteins'], gtf, fasta ] }
        .multiMap { meta, gtf_file, fasta_file ->
            gtf_input: [ meta, gtf_file ]
            fasta_input: fasta_file
        }

    PROTEIN_EXTRACT (
        ch_gtf_fasta_inputs_prot.gtf_input,
        ch_gtf_fasta_inputs_prot.fasta_input
    )

    BLAST_MAKEBLASTDB ( PROTEIN_EXTRACT.out.gffread_fasta.map { meta, fa -> [ meta, fa ] } )
    ch_blast_protein_db = BLAST_MAKEBLASTDB.out.db
//    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)


    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    fai              = ch_fai                    // channel: path(genome.fai)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    rrna_fastas      = ch_rrna_fastas            // channel: path(sortmerna_fasta_list)
    sortmerna_index  = ch_sortmerna_index        // channel: path(sortmerna/index/)
    star_index       = ch_star_index             // channel: path(star/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
    cds_fasta        = ch_cds_fasta
    mrna_fasta       = ch_mrna_fasta
    lncrna_fasta     = ch_lncrna_fasta
    known_lncrna_gtf = ch_known_lncrna_gtf
    ncrna_fasta      = ch_ncrna_fasta
    blast_protein_db = ch_blast_protein_db
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
