process FASTQC {
	publishDir "$params.outdir/1_preprocessed_reads/1_quality_check",	pattern: "*.html",	mode: "copy", overwrite: false
	publishDir "$params.outdir/1_preprocessed_reads/1_quality_check",	pattern: "*.zip",	mode: "copy", overwrite: false
    input:
        path reads
    output:
        path "*.html",	emit: fastqc_html
        path "*.zip",	emit: fastqc_zip
    script:
        """
        #!/bin/bash
        fastqc \
        -o \$PWD \
        --noextract \
        ${reads}
        """
}

process TRIMGALORE {
    tag "$sample_id"
    publishDir "$params.outdir/1_preprocessed_reads/2_trimmed_reads", pattern: "trimmed_output/*.html",   mode: 'copy',   overwrite: false
    publishDir "$params.outdir/1_preprocessed_reads/2_trimmed_reads", pattern: "trimmed_output/.zip",     mode: 'copy',   overwrite: false
    publishDir "$params.outdir/1_preprocessed_reads/2_trimmed_reads", pattern: "trimmed_output/.txt",     mode: 'copy',   overwrite: false
    publishDir "$params.outdir/1_preprocessed_reads/2_trimmed_reads", pattern: "trimmed_output/*.fq.gz",  mode: 'copy',   overwrite: false
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("trimmed_output/*val_1.fq.gz")      , emit: trimmed_read1
        tuple val(sample_id), path("trimmed_output/*val_2.fq.gz")      , emit: trimmed_read2
        path "trimmed_output/*.html"            , emit: html
        path "trimmed_output/*.zip"             , emit: zip
        path "trimmed_output/*report.txt"       , emit: log
    script:
        """
        #!/bin/bash
        mkdir trimmed_output
        trim_galore \
        --paired \
        --fastqc \
        --quality 20 \
        --length 25 \
        --cores ${task.cpus} \
	--gzip \
        -o trimmed_output \
        ${reads[0]} ${reads[1]}
        """
}

process REVERSECOMPLEMENT {
    tag "$sample_id"
    publishDir "$params.outdir/1_preprocessed_reads/3_revComp", mode: "copy"

    input:
        tuple val(sample_id), path(read1), path(read2)
    output:
        tuple val(sample_id), path("${sample_id}_R1_final.fq.gz"), path("${sample_id}_R2_final.fq.gz"), emit: reversed_reads

    script:
    """
    #!/bin/bash
    if [[ "$params.lib_type" == "rf" ]]; then
        seqkit seq -r -p ${read1} | gzip > ${sample_id}_R1_final.fq.gz
        ln -s ${read2} ${sample_id}_R2_final.fq.gz
    elif [[ "$params.lib_type" == "fr" ]]; then
        ln -s ${read1} ${sample_id}_R1_final.fq.gz
        seqkit seq -r -p ${read2} | gzip > \${sample_id}_R2_final.fq.gz
    else
        ln -s ${read1} ${sample_id}_R1_final.fq.gz
        ln -s ${read2} ${sample_id}_R2_final.fq.gz
    fi
    """
}

process MINRECURRENCE {
    input: 
        path metafile
    output:
        stdout emit: min_recurrence
    script:
    """
    #!/bin/bash
    awk '{print \$1}' $metafile | sort | uniq -c | sort -nr | head -n 1 | awk '{print \$1}' | awk '{printf "%s", \$0}'
    """
}

process JFCOUNT {
    publishDir "$params.outdir/2_kmer_processing/1_jellyfish/count",  pattern: "*.jf",    mode: "copy",   overwrite: false

    tag { "${condition}_${replicate}" }

    input:
    tuple val(sample_id), val(condition), val(replicate), path(read1), path(read2)
    
    output:
        path("${condition}_${replicate}.jf"),    emit: jellyfish_count
    script:
        """
        #!/bin/bash
        jellyfish count -m 31 -s 100000000 -t ${task.cpus} -o ${condition}_${replicate}.jf -F 2 <(zcat ${read1}) <(zcat ${read2})
        """
}

process JFDUMP {
    publishDir "$params.outdir/2_kmer_processing/1_jellyfish/dump",  pattern: "*.txt.gz",    mode: "copy",   overwrite: false
    input:
        path jf_count
    output:
        path "*.txt.gz",   emit: jellyfish_dump
    script:
        """
        #!/bin/bash
	set -euo pipefail
	TMP=\$PWD/tmp
        mkdir -p \$TMP
        jellyfish dump -c $jf_count \\
	| sort --parallel=${task.cpus} -S ${task.memory.toGiga()-2}G -T \$TMP -k1 \\
	| gzip -c > ${jf_count.simpleName}.txt.gz 
        """
}

process JOINCOUNTS {
	publishDir "$params.outdir/2_kmer_processing/1_jellyfish/join", pattern: "*.tsv.gz", mode: "copy",   overwrite: true
	input:
        path jf_output
        val min_recurrence
	output: 
	    path "*.tsv.gz", emit: joined_kmer
	script:
	    """
        #!/bin/bash
        echo -n "tag" | gzip -c > kmer-counts.tsv.gz
        for s in ${jf_output}
        do
	        echo -ne "\t"\$s | gzip -c >> kmer-counts.tsv.gz
        done
        echo "" | gzip -c >> kmer-counts.tsv.gz
        joinCounts -r $min_recurrence -a 5 $jf_output | gzip -c >> kmer-counts.tsv.gz
	    """
}

process KAMRATINDEX {
    publishDir "$params.outdir/2_kmer_processing/2_kamrat/kamrat_index",  mode: "copy",   overwrite: false
    input:
        path joincounts_output
    output:
        path "*",   emit: kamrat_index
    script:
        """
        #!/bin/bash
        mkdir kamrat_idx
        kamrat index -intab $joincounts_output -outdir kamrat_idx -klen 31
        """
}

process KAMRATMASK {
    publishDir "$params.outdir/2_kmer_processing/2_kamrat",  mode: "copy",  overwrite: true
    input:
        path kamrat_index
	path transcriptome
    output:
        path "*.bin",    emit: masked_kmers
        path "*.log",    emit: masking_log
    script:
        """
        #!/bin/bash
        kamrat mask -idxdir $kamrat_index -fasta $params.transcriptome \
        -outpath masked_kmers.bin &> masking.log
        """
}

process KAMRATMERGE {
    publishDir "$params.outdir/2_kmer_processing/2_kamrat", mode: "copy",   overwrite: true
    input:
        path kamrat_index
        path masked_kmers
    output:
        path "*.tsv",   emit: masked_merged_kmers
        path "*.log",   emit: merge_log
    script:
        """
        #!/bin/bash
        kamrat merge -idxdir $kamrat_index \
        -overlap 30-15 -interv none -with $masked_kmers -withcounts mean -outpath masked_merged_kmers.tsv &> merge.log
        """
}

process ANNOTCONFIG_FIRST {
    input:
        path contigs_for_annotation
        path genome
        path annotation
        path repeat
        val(lib_type)
    output:
        path "*.json",  emit: annot_config_first
    script:	
        """
        #!/bin/bash
        
        echo "{
            \\"mode\\": \\"index\\",
	        \\"preset\\": \\"splice\\",
            \\"input_file\\": \\"$contigs_for_annotation\\",
            \\"sequence_col\\": \\"contig\\",
            \\"id_col\\": \\"tag\\",
            \\"map_to\\" :[\\"human\\"],
            \\"supp_map_to\\":[\\"repeats\\"],
            \\"supp_map_to_fasta\\" :[\\"$repeat\\"],
            \\"output_dir\\": \\"output\\",
            \\"keep_col\\": \\"all\\",
            \\"library_type\\": \\"$lib_type\\",
            \\"reference\\": \\"$genome\\",
            \\"annotation\\": \\"$annotation\\",
            \\"threads\\": 8,
            \\"star_index\\": \\"\\",
            \\"minimap2_index\\" : \\"\\"
        }" >> config.json
        """
}


process ANNOTCONTIGS_FIRST {
    publishDir "$params.outdir/3_contig_annotation", mode: "copy", overwrite: true
    input:
        path contigs_for_annotation
        path annot_config_first
        path genome
        path annotation
        path repeat
	path snakefile
	path script
    output:
        path "output",                          emit: annotation_outputs_first
        path "output/merged_annotation.tsv",    emit: annotated_contigs_first
        path "human_index/star",		emit: star_index
        path "human_index/minimap2",	        emit: minimap2_index
    script:
        """
        #!/bin/bash
        export PATH=/opt/conda/envs/AnnotateContig/bin:$PATH
        snakemake -s $snakefile --configfile $annot_config_first --cores ${task.cpus}
        """
}

process ANNOTCONFIG {
    input:
        path contigs_for_annotation
        path genome
        path annotation
        path repeat
        val(lib_type)
        path human_index
    output:
        path "*.json",  emit: annot_config
    script:	
        """
        #!/bin/bash
            echo "{
                \\"mode\\": \\"table\\",
		        \\"preset\\": \\"splice\\",
                \\"input_file\\": \\"$contigs_for_annotation\\",
                \\"sequence_col\\": \\"contig\\",
                \\"id_col\\": \\"tag\\",
                \\"map_to\\" :[\\"human\\"],
                \\"supp_map_to\\":[\\"repeats\\"],
                \\"supp_map_to_fasta\\" :[\\"$repeat\\"],
                \\"output_dir\\": \\"output\\",
                \\"keep_col\\": \\"all\\",
                \\"library_type\\": \\"$lib_type\\",
                \\"reference\\": \\"$genome\\",
                \\"annotation\\": \\"$annotation\\",
		        \\"threads\\": 8,
                \\"star_index\\": \\"$human_index/star\\",
                \\"minimap2_index\\" : \\"$human_index/minimap2/human.mmi\\"
            }" >> config.json
        """
}

process ANNOTCONTIGS {
    publishDir "$params.outdir/3_contig_annotation", mode: "copy", overwrite: true
    input:
        path contigs_for_annotation
        path annot_config
        path genome
        path annotation
        path repeat
	    path snakefile
	    path script
        path human_index
    output:
        path "output",                          emit: annotation_outputs
        path "output/merged_annotation.tsv",    emit: annotated_contigs
    script:
        """
        #!/bin/bash
        export PATH=/opt/conda/envs/AnnotateContig/bin:$PATH
        snakemake -s $snakefile --configfile $annot_config --cores ${task.cpus}
        """
}

process STAR{
    publishDir "$params.outdir/4_aligned_sequences/bams",   pattern: "*sortedByCoord.out.bam",  mode: 'copy', overwrite: false
    publishDir "$params.outdir/4_aligned_sequences/logs",   pattern: "*Log*",                   mode: 'copy', overwrite: false
    publishDir "$params.outdir/4_aligned_sequences/logs",   pattern: "*out.tab*",               mode: 'copy', overwrite: false
    input:
        tuple val(sample_id), path(trimmed_reads1), path(trimmed_reads2)
        path star_index
    output:
        path "*sortedByCoord.out.bam"   , emit: sorted_bam
        path "*Log.final.out"           , emit: log_final
        path "*Log.out"                 , emit: log_out
        path "*Log.progress.out"        , emit: log_progress
        path "*ReadsPerGene.out.tab"    , emit: readspergene_out
        path "*SJ.out.tab"              , emit: splice_junctions
    script:
        def prefix = trimmed_reads1.name.replaceAll(/_val_1\.f(ast)?q(\.gz)?$/, '')
        """
        #!/bin/bash
        STAR \
        --genomeDir $star_index \
        --readFilesCommand gunzip -c \
        --readFilesIn $trimmed_reads1 $trimmed_reads2 \
        --outFileNamePrefix ${prefix}_ \
        --outSAMtype BAM SortedByCoordinate \
	    --quantMode GeneCounts \
        --runThreadN ${task.cpus} \
        --outSJfilterOverhangMin 12 12 12 12 \
        --outSJfilterCountUniqueMin 1 1 1 1 \
        --outSJfilterCountTotalMin 1 1 1 1 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --scoreGapNoncan 0 \
        --scoreGapGCAG 0 \
        --scoreGapATAC 0 
        """
}

process EXTRACTJUNCTIONS{
    input:
        path splice_junctions
    output:
        path "*_junctions.bed"  , emit: merged_splices
    script:
        """
        #!/bin/bash
        awk -F'\t' '{
            for (i=2; i<=NF; i++) {
                if (\$i ~ /^[0-9]+\$/) {
                    \$i = \$i+0
                }
            }       
            if (\$6 != 1 && \$7 >= 3 || \$8 >= 5 && \$4 != 0) {
                if (\$4 == 1) { \$4 = "+" }
                else if (\$4 == 2) { \$4 = "-" }
                print \$1"\t"\$2-1"\t"\$3"\t"\$4
            }
        }' $splice_junctions > ${splice_junctions.simpleName}_junctions.bed
        """ 
}


process MERGESPLICES{
    publishDir "$params.outdir/5_candidates", pattern: "*final_merged_splices.txt",   mode: 'copy',   overwrite: true
    input:
        path merged_splices
    output:
        path "*final_merged_splices.txt",   emit: final_splices
    script:
        """
        #!/bin/bash
        echo "$merged_splices" | tr -d '[]' | tr ', ' '\n' | while read -r file; do
            cat "\$file"
        done | sort -k1,1V -k2,2n -k3,3n | uniq > final_merged_splices.txt
        """
}

process INTRONRETENTIONS {
    publishDir "$params.outdir/5_candidates", pattern: "retained_introns.bed",     mode: 'copy', overwrite: true
    publishDir "$params.outdir/5_candidates", pattern: "intron_coverage_full.tsv", mode: 'copy', overwrite: true
    cpus 8
    memory '4 GB'
    input:
        path intron_file
        path sorted_bam
    output:
        path "retained_introns.bed",      emit: retained_introns
        path "intron_coverage_full.tsv",  emit: intron_coverage
    script:
    """
    #!/bin/bash
    set -euo pipefail
    THREADS=${task.cpus}

    # merge all bams
    samtools merge -@ \$THREADS -f merged.bam ${sorted_bam}
    samtools index -@ \$THREADS merged.bam

    # split reads by strand
    build_strand_bam () {
        samtools view -@ \$THREADS -b \$2 merged.bam > _p1.bam
        samtools view -@ \$THREADS -b \$3 merged.bam > _p2.bam
        samtools merge -@ \$THREADS -f \$1 _p1.bam _p2.bam
        samtools index -@ \$THREADS \$1
        rm -f _p1.bam _p2.bam
    }

    if [ "$params.lib_type" == "rf" ]; then
        build_strand_bam plus.bam  "-f 128 -F 16" "-f 80"
        build_strand_bam minus.bam "-f 144"       "-f 64 -F 16"
    elif [ "$params.lib_type" == "fr" ]; then
        build_strand_bam plus.bam  "-f 64 -F 16"  "-f 144"
        build_strand_bam minus.bam "-f 80"        "-f 128 -F 16"
    else
        echo "ERROR: lib_type '$params.lib_type' desteklenmiyor (rf|fr)" >&2
        exit 1
    fi

    # split introns by strand
    awk -F'\\t' '\$6=="+"' ${intron_file} > introns.plus.bed
    awk -F'\\t' '\$6=="-"' ${intron_file} > introns.minus.bed

    # calculate strand-mathed breadth 
    mosdepth -t \$THREADS -n --by introns.plus.bed  --thresholds 1 plus_md  plus.bam
    mosdepth -t \$THREADS -n --by introns.minus.bed --thresholds 1 minus_md minus.bam

    #    thresholds.bed.gz columns: #chrom start end region 1X
    {
      echo -e "chr\\tstart\\tend\\tintron_id\\tcovered_bp\\tlength\\tfraction\\tstrand"
      zcat plus_md.thresholds.bed.gz  | awk 'BEGIN{OFS="\\t"} \$1!~/^#/ {len=\$3-\$2; f=(len>0?\$5/len:0); print \$1,\$2,\$3,\$4,\$5,len,f,"+"}'
      zcat minus_md.thresholds.bed.gz | awk 'BEGIN{OFS="\\t"} \$1!~/^#/ {len=\$3-\$2; f=(len>0?\$5/len:0); print \$1,\$2,\$3,\$4,\$5,len,f,"-"}'
    } > intron_coverage_full.tsv

    # retained = >%85 of introns >=1x covered
    awk -F'\\t' 'NR>1 && \$7>0.85 {print \$4}' intron_coverage_full.tsv | sort -u > retained_introns.bed
    """
}

process UNMAPPEDCONTIGS {
    publishDir "$params.outdir/5_candidates", pattern: "*not_covered_tags.txt", mode: 'copy', overwrite: true

    input:
        path annotated_contigs
    output:
        path "merged_annotation_filtered.tsv",  emit: filtered_contigs
    script:
        """
        # Filter annotated_contigs to skip rows with 'mapped_to' == "NA"
        awk -F'\\t' '
        NR==1 {
            # Find which column is mapped_to
            for (i=1; i<=NF; i++) {
                if(\$i=="mapped_to") col=i
            }
            # Print the header line and move to next line
            print \$0
            next
        }
        # Print only lines where the mapped_to column is not "NA"
        \$col != "NA" {
            print \$0
        }
        ' $annotated_contigs > merged_annotation_filtered.tsv
        """
}
	

process EXONEXTENSIONS{
    publishDir "$params.outdir/5_candidates", pattern: "*extended_exons.bed",   mode: 'copy',   overwrite: true
    input:
        path input_bed
        path sorted_bam
    output:
        path "*extended_exons.bed",   emit: extended_exons
    script:
        // Decide which AWK command to use based on existence of exon file
        def awkCommand
        if (params.exon_file) {
            awkCommand = """
                awk '{start1=\$2-10; end1=\$2; start2=\$3; end2=\$3+10; print \$1 \"\\t\" start1 \"\\t\" end1 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6; print \$1 \"\\t\" start2 \"\\t\" end2 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6;}' ${input_bed} | sort -k1,1V -k2,2n -k3,3n > output.bed
            """
        } else {
            awkCommand = """
                awk '{start1=\$2; end1=start1+10; start2=\$3-10; end2=\$3; print \$1 \"\\t\" start1 \"\\t\" end1 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6; print \$1 \"\\t\" start2 \"\\t\" end2 \"\\t\" \$4 "\\t" \$5 \"\\t\" \$6;}' ${input_bed} | sort -k1,1V -k2,2n -k3,3n > output.bed
            """
        }
        """
        #!/bin/bash
        bam_list=($sorted_bam)

        # awk '{start1=\$2; end1=start1+10; start2=\$3-10; end2=\$3; print \$1 \"\\t\" start1 \"\\t\" end1 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6; print \$1 \"\\t\" start2 \"\\t\" end2 \"\\t\" \$4 \"\\t\" \$5 \"\\t\" \$6;}' ${input_bed} | sort -k1,1V -k2,2n -k3,3n > output.bed
        # Execute the chosen AWK command
        ${awkCommand}

        bedtools_command="bedtools coverage -a output.bed"

        for bam in \${bam_list[@]}; do
            bedtools_command+=" -b \$bam"
        done

        bedtools_command+=" -s -sorted -split -counts \\
        | awk 'BEGIN {OFS=FS=\\"\\t\\"}; {if (\\\$7 > 5) print}' \\
        | cut -f 1-6 | uniq > extended_exons.bed"

        echo "Executing command:"
        echo "\$bedtools_command"
        eval "\$bedtools_command"
        """
}

process BAMINDEX {
    publishDir "$params.outdir/4_aligned_sequences/bams", mode: 'copy', overwrite: false
    input:
        path sorted_bam
    output:
        path "*.bai"      ,emit: sorted_bam_bai
    script:
        """
        #!/bin/bash
        samtools index $sorted_bam
        """
}


process FINDCANDIDATES_DE {
    publishDir "$params.outdir/5_candidates", mode: "copy", overwrite: true
    input:
        path cand_finder
        path filtered_contigs
        path intron_file
        path metafile
        val(control)
        val(test)
        path final_splices
        path retained_introns
        path extended_exons
        path gtf
        path functions
        path annotation_outputs
    output:
        path "*", emit: candidate_results, optional: true
        path "*_candidates.tsv",    emit: filtered_introns
        path "conditions.Rdata",    emit: conditions
    script:
        """
        #!/bin/bash
        Rscript $cand_finder $filtered_contigs $intron_file $metafile $control $test $final_splices $retained_introns $extended_exons $gtf
        """
}

process FINDCANDIDATES {
    publishDir "$params.outdir/5_candidates", mode: "copy", overwrite: true
    input:
        path cand_finder
        path filtered_contigs
        path intron_file
        path metafile
        val(sample_name)
        path final_splices
        path retained_introns
        path extended_exons
        path gtf
        path functions
        path annotation_outputs
    output:
        path "*", emit: candidate_results, optional: true
        path "*_candidates.tsv",    emit: filtered_introns
        path "conditions.Rdata",    emit: conditions
    script:
        """
        #!/bin/bash
        Rscript $cand_finder $filtered_contigs $intron_file $metafile $sample_name $final_splices $retained_introns $extended_exons $gtf
        """
}

process COVEREDCANDIDATES {
    publishDir "$params.outdir/5_candidates", pattern: "not_covered_introns.txt", mode: 'copy', overwrite: true
    input:
        path filtered_introns
        path sorted_bam
    output:
        path "not_covered_introns.txt",    emit: not_covered_introns
        path "filtered_introns.bed",    emit: filtered_introns_bed
    script:
        """
        #!/bin/bash
        awk -F'\\t' '
            NR==1 {
                # Identify the column positions of the headers
                for (i = 1; i <= NF; i++) {
                    if (\$i == "chr") col_chr = i
                    else if (\$i == "start") col_start = i
                    else if (\$i == "end") col_end = i
                    else if (\$i == "intron_id") col_intron_id = i
                    else if (\$i == "strand") col_strand = i
                }
                # Skip printing the header row itself
                next
            }
            # Print only rows where the chromosome column has "chr"
            (\$col_chr ~ /chr/) {
                # Print the BED columns, inserting a 0 before strand
                print \$col_chr "\\t" \$col_start "\\t" \$col_end "\\t" \$col_intron_id "\\t0\\t" \$col_strand
            }
        ' OFS='\\t' $filtered_introns | sort -k1,1V -k2,2n -k3,3n > filtered_introns.bed

        # Build the bedtools coverage command
        bam_list=($sorted_bam)
        bedtools_command="bedtools coverage -a filtered_introns.bed"

        for bam in \${bam_list[@]}; do
            bedtools_command+=" -b \$bam"
        done

        # 4) Execute coverage and capture tags with coverage < 5
        bedtools_command+=" -s -sorted -split -counts |
        awk 'BEGIN {OFS=FS=\\"\\t\\"}; {if (\\\$7 < 5) print}' |
        cut -f4 | uniq > not_covered_introns.txt"

        echo "Executing command:"
        echo "\$bedtools_command"
        eval "\$bedtools_command"
        """
}

process FINALFILTERING {
    publishDir "$params.outdir/5_candidates", mode: "copy", overwrite: true
    input:
        path coverage_filter
        path filtered_introns
        path not_covered_introns
        path conditions
    output: 
        path "*_covered_candidates.tsv",    emit: final_intron_candidates
    script:
        """
        #!/bin/bash 
        Rscript $coverage_filter $filtered_introns $not_covered_introns
        """

}
    


