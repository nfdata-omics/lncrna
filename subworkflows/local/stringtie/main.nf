include { STRINGTIE_STRINGTIE }         from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE }             from '../../../modules/nf-core/stringtie/merge/main'
include { GFFCOMPARE }                  from '../../../modules/nf-core/gffcompare/main'
include { FILTER_GTF_BY_CLASSCODE }     from '../../../modules/local/filter_gtf_by_classcode'
include { GTF_TO_FASTA }                from '../../../modules/local/gtf_to_fasta'


workflow STRINGTIE_WORKFLOW {
    take:
        bam_sorted
        ch_gtf
        ch_fasta
        ch_fai

    main:
        ch_versions = Channel.empty()
        ch_stringtie_gtf = Channel.empty()
        //
        // STRINGTIE: Transcript assembly and quantification for RNA-SeQ
        //
        STRINGTIE_STRINGTIE(bam_sorted, ch_gtf)
        ch_stringtie_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf.map { meta, transcript_gtf -> [ transcript_gtf ] }.collect()
//        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

        //
        // STRINGTIE: Merge transcript assemblies into a non-redundant annotation
        //
        STRINGTIE_MERGE (ch_stringtie_gtf, ch_gtf)
        ch_stringtie_gtf_merged = STRINGTIE_MERGE.out.gtf
            .map { gtf_file ->
                def meta = [ id: 'stringtie_merged' ]
                return [ meta, gtf_file ]
            }
//        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)


        //
        // GFFCOMPARE: Annotate with class codes
        //
        // Prepare FASTA with FAI index
        ch_fasta_meta_fai = ch_fasta.combine(ch_fai)
            .map { fasta, fai ->
                def meta2 = [ id: fasta.getBaseName(), description: 'Genome FASTA with index' ]
                return [ meta2, fasta, fai ]
            }

        // Prepare reference GTF
        ch_reference_gtf = ch_gtf.map { gtf_file ->
            def meta3 = [ id: gtf_file.getBaseName(), description: 'Reference GTF file' ]
            return [ meta3, gtf_file ]
        }

        GFFCOMPARE (
            ch_stringtie_gtf_merged,
            ch_fasta_meta_fai,
            ch_reference_gtf
        )
//        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

        //
        // FILTER: lncRNA candidates by class code (i, u, x)
        //
        FILTER_GTF_BY_CLASSCODE ( GFFCOMPARE.out.annotated_gtf.map { it[1] } )
        ch_lncrna_candidates = FILTER_GTF_BY_CLASSCODE.out.filtered_gtf
            .map { gtf_file ->
                def meta = [ id: 'lncrna_candidates' ]
                return [ meta, gtf_file ]
            }
        ch_classcode_stats = FILTER_GTF_BY_CLASSCODE.out.stats
//        ch_versions = ch_versions.mix(FILTER_GTF_BY_CLASSCODE.out.versions)

        //
        // CONVERT: GTF to FASTA
        //
        GTF_TO_FASTA ( ch_lncrna_candidates, ch_fasta )
        ch_lncrna_fasta = GTF_TO_FASTA.out.fasta
//        ch_versions = ch_versions.mix(GTF_TO_FASTA.out.versions)


    emit:
        stringtie_gtf_merged = ch_stringtie_gtf_merged
        lncrna_candidates    = ch_lncrna_candidates                       // Filtered GTF by class codes
        classcode_stats      = ch_classcode_stats
        lncrna_fasta         = ch_lncrna_fasta                            // FASTA sequences
        gffcompare_annotated = GFFCOMPARE.out.annotated_gtf               // GTF with class codes
        tmap                 = GFFCOMPARE.out.tmap.map { meta, tmap -> tmap }
        versions             = ch_versions
}
