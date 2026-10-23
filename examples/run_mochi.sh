#!/bin/bash
# MOCHI run script. Edit the paths below, then run: bash run_mochi.sh

SIF=/path/to/mochi.sif
DATA=/path/to/data                  # parent directory of all inputs (bound into the container)
REF=$DATA/references

# Mode: candidate discovery only (single condition in the metafile).

# For differential expression, use metafile_DE.txt and set:
# MODE="--DE --control_group CTRL --test_group TEST"

MODE="--sample_name CTRL"

# Indexes: STAR and minimap2 needs index files. --index builds them on the first run. On later runs, reuse them:
# INDEX="--human_index /path/to/human_index --star_index /path/to/human_index/star"
INDEX="--index"

apptainer run --cleanenv -B "$DATA" "$SIF" \
    -resume \
    --metafile      metafile_single.txt \
    --raw_reads     "$DATA/fastq/*{_R1,_R2}.fastq.gz" \
    --genome        "$REF/GRCh38.p14.genome.fa.gz" \
    --transcriptome "$REF/gencode.v47.transcripts.fa" \
    --annotation    "$REF/gencode.v47.annotation.gtf.gz" \
    --exon_file     "$REF/gencode_v47_exons.bed" \
    --intron_file   "$REF/gencode_v47_introns.bed" \
    --repeat        "$REF/human_repeat_ref.fasta" \
    --lib_type      rf \
    $INDEX \
    $MODE
