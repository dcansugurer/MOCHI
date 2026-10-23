// Input parameters from CLI

params.raw_reads 	    = ""
params.trimmed_reads 	= ""
params.metafile 	    = ""
params.transcriptome 	= ""
params.genome 		    = ""
params.annotation 	    = ""
params.repeat 		    = ""
params.human_index 	    = ""
params.star_index 	    = ""
params.minimap2_index 	= ""
params.intron_file 	    = ""
params.lib_type 	    = ""
params.sample_name 	    = ""
params.control_group 	= ""
params.test_group 	    = ""
params.index 		    = false
params.joined_counts 	= false
params.exon_file 	    = false
params.DE 		        = false

// Default Parameters

params.outdir       = "${launchDir}/output"
cand_finder_with_DE = Channel.of("$projectDir/script/cand_finder_with_DE.R")
cand_finder_no_DE   = Channel.of("$projectDir/script/cand_finder_no_DE.R")
functions           = Channel.of("$projectDir/script/functions.R")
coverage_filter     = Channel.of("$projectDir/script/filter_covered.R")

// Channels
Channel.fromPath(params.transcriptome)
.set { transcriptome_ch}

Channel.fromPath(params.genome)
.set { genome_ch }

Channel.fromPath(params.annotation)
.set { annotation_ch }

Channel.fromPath(params.repeat)
.set { repeat_ch }

Channel.of(params.lib_type)
.set { lib_type }

Channel.fromPath(params.intron_file)
.set { intron_file }

log.info """\
                         ...''''''''...
                     ..''              ''..
                   .'                      '.
                 .'                          '.
                .                              .
               #                                #
              #                                  #
             .                                    .
            .       .'#.            .'#.           .
            #       ####  .      .  ####           #
           .     . .       .    .      . .          .
           '    '.'.'.      ....      '.'.'.        '.
          #       ' '                   ' '          #
          #                                         '#
          #.                                       .'#
           #'.                                   .'.#
            '.'.'.                           . .'.#'
              '.#.'.'.'. . . . . . . . .'.'.'.'.#'
                 ''.#.'.'.'.'.'.'.'.'.'.'.#.'''
                       '''''''''''''''''

          ##     ##  #######   ######  ##     ## ####
          ###   ### ##     ## ##    ## ##     ##  ##
          #### #### ##     ## ##       ##     ##  ##
          ## ### ## ##     ## ##       #########  ##
          ##     ## ##     ## ##       ##     ##  ##
          ##     ## ##     ## ##    ## ##     ##  ##
          ##     ##  #######   ######  ##     ## ####

         Mapping of Cryptic Hidden Intronic Transcripts

         reads        : ${params.raw_reads}
         transcriptome: ${params.transcriptome}
         genome       : ${params.genome}
         annotation   : ${params.annotation}
         repeat       : ${params.repeat}
         outdir       : ${params.outdir}
         """
         .stripIndent()

nextflow.enable.dsl=2

include { FASTQC }              from './modules.nf'
include { TRIMGALORE }          from './modules.nf'
include { REVERSECOMPLEMENT }   from './modules.nf'
include { MINRECURRENCE }       from './modules.nf'
include { JFCOUNT }             from './modules.nf'
include { JFDUMP }              from './modules.nf'
include { JOINCOUNTS }          from './modules.nf'
include { KAMRATINDEX }         from './modules.nf'
include { KAMRATMASK }          from './modules.nf'
include { KAMRATMERGE }         from './modules.nf'
include { ANNOTCONFIG_FIRST }   from './modules.nf'
include { ANNOTCONTIGS_FIRST }  from './modules.nf'
include { ANNOTCONFIG }         from './modules.nf'
include { ANNOTCONTIGS }        from './modules.nf'
include { STAR }                from './modules.nf'
include { EXTRACTJUNCTIONS }    from './modules.nf'
include { MERGESPLICES }        from './modules.nf'
include { INTRONRETENTIONS }    from './modules.nf'
include { UNMAPPEDCONTIGS }     from './modules.nf'
include { EXONEXTENSIONS }      from './modules.nf'
include { BAMINDEX }            from './modules.nf'
include { FINDCANDIDATES }      from './modules.nf'
include { FINDCANDIDATES_DE }   from './modules.nf'
include { COVEREDCANDIDATES }   from './modules.nf'
include { FINALFILTERING }      from './modules.nf'




workflow {
    Channel.fromPath(params.metafile)
        | splitCsv(header: true, sep: '\t')
        | map { row -> tuple(row.sample_id, row) }
        | set { meta_ch }
    
    MINRECURRENCE(Channel.fromPath(params.metafile))
    
    if ( params.joined_counts ) {
        KAMRATINDEX(Channel.fromPath(params.joined_counts))
        
        Channel.fromFilePairs(params.trimmed_reads)
        .set { trimmed_reads_ch }
        
    } else {
        Channel.fromFilePairs(params.raw_reads)
        .set { read_pairs_ch } 

        Channel.fromPath(params.raw_reads)
        .set { reads_ch }
        FASTQC(reads_ch)
        TRIMGALORE(read_pairs_ch)
        TRIMGALORE
            .out
            .trimmed_read1
            .combine(TRIMGALORE.out.trimmed_read2, by:0)
            .set { trimmed_reads_ch }
        REVERSECOMPLEMENT(TRIMGALORE.out.trimmed_read1.combine(TRIMGALORE.out.trimmed_read2, by:0))

        REVERSECOMPLEMENT.out.map { [it[0], [it[1], it[2]]] }
        .set { revComp_files }
        
        meta_ch
        .join(revComp_files)
        .map { sample_id, meta, reads -> [sample_id, meta.condition, meta.replicate, reads[0], reads[1]] }
        .set { samples }


        JFCOUNT(samples)
        JFDUMP(JFCOUNT.out.jellyfish_count)
        JOINCOUNTS(JFDUMP.out.jellyfish_dump.collect(), MINRECURRENCE.out.min_recurrence)
        KAMRATINDEX(JOINCOUNTS.out.joined_kmer)
    }

    KAMRATMASK(KAMRATINDEX.out.kamrat_index, Channel.fromPath(params.transcriptome))
    KAMRATMERGE(KAMRATINDEX.out.kamrat_index, KAMRATMASK.out.masked_kmers)

    Channel.of(file("$projectDir/annotate_contigs/Snakefile")).set { snakefile_ch }
    Channel.of(file("$projectDir/annotate_contigs/script")).set { script_ch }

    if ( params.index ) {
        ANNOTCONFIG_FIRST(KAMRATMERGE.out.masked_merged_kmers, genome_ch, annotation_ch, repeat_ch, lib_type)
        ANNOTCONTIGS_FIRST(KAMRATMERGE.out.masked_merged_kmers, ANNOTCONFIG_FIRST.out.annot_config_first, 
        genome_ch, annotation_ch, repeat_ch, snakefile_ch, script_ch)
        ANNOTCONTIGS_FIRST
            .out
            .annotated_contigs_first
            .set { annotated_contigs_ch }
        ANNOTCONTIGS_FIRST
            .out
            .annotation_outputs_first
            .set { annotation_outputs_ch }
        ANNOTCONTIGS_FIRST
            .out
            .star_index
            .set { star_index_ch }
        ANNOTCONTIGS_FIRST
            .out
            .minimap2_index
            .set { minimap2_index_ch }
    }
    else {
        Channel.fromPath(params.human_index)
            .set { human_index_ch }
        ANNOTCONFIG(KAMRATMERGE.out.masked_merged_kmers, genome_ch, annotation_ch, repeat_ch, lib_type, human_index_ch)
        ANNOTCONTIGS(KAMRATMERGE.out.masked_merged_kmers, ANNOTCONFIG.out.annot_config, 
        genome_ch, annotation_ch, repeat_ch, snakefile_ch, script_ch, human_index_ch)
        ANNOTCONTIGS
            .out
            .annotated_contigs
            .set { annotated_contigs_ch }
        ANNOTCONTIGS
            .out
            .annotation_outputs
            .set { annotation_outputs_ch }
    }

    if (params.star_index) {
        // If STAR index is provided, use it directly
        Channel.fromPath(params.star_index)
            .set { star_index_ch }
    } else {
       ANNOTCONTIGS_FIRST
            .out
            .star_index
            .set { star_index_ch }
    }

    if (params.exon_file) {
        // If user provides exon file, calculate more accurate exon extensions
        Channel.fromPath(params.exon_file)
            .set { input_bed_ext }
    } else {
        Channel.fromPath(params.intron_file)
            .set { input_bed_ext }
    }

    STAR(trimmed_reads_ch, star_index_ch.first())
    EXTRACTJUNCTIONS(STAR.out.splice_junctions)
    MERGESPLICES(EXTRACTJUNCTIONS.out.merged_splices.collect())
    INTRONRETENTIONS(intron_file, STAR.out.sorted_bam.collect())
    UNMAPPEDCONTIGS(annotated_contigs_ch)
    EXONEXTENSIONS(input_bed_ext, STAR.out.sorted_bam.collect())
    BAMINDEX(STAR.out.sorted_bam)

    if (params.DE) {
    FINDCANDIDATES_DE(cand_finder_with_DE, UNMAPPEDCONTIGS.out.filtered_contigs,
        intron_file, Channel.fromPath(params.metafile), Channel.of(params.control_group), Channel.of(params.test_group),
        MERGESPLICES.out.final_splices, INTRONRETENTIONS.out.retained_introns, EXONEXTENSIONS.out.extended_exons, 
        annotation_ch, functions, annotation_outputs_ch)
    filtered_introns = FINDCANDIDATES_DE.out.filtered_introns
    conditions       = FINDCANDIDATES_DE.out.conditions
    } else {
    FINDCANDIDATES(cand_finder_no_DE, UNMAPPEDCONTIGS.out.filtered_contigs,
        intron_file, Channel.fromPath(params.metafile), Channel.of(params.sample_name),
        MERGESPLICES.out.final_splices, INTRONRETENTIONS.out.retained_introns, EXONEXTENSIONS.out.extended_exons, 
        annotation_ch, functions, annotation_outputs_ch)
    filtered_introns = FINDCANDIDATES.out.filtered_introns
    conditions       = FINDCANDIDATES.out.conditions
    }
    COVEREDCANDIDATES(filtered_introns, STAR.out.sorted_bam.collect())
    FINALFILTERING(coverage_filter, filtered_introns, COVEREDCANDIDATES.out.not_covered_introns, conditions)
}
