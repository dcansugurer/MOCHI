#------------------------------------ Required Functions for MOCHI Script ------------------------------------

# Function for reading contig table
read_contig_file = function(contig_file) {
  # Read tab separated contig annotation
  contigs = read.table(contig_file, header = TRUE, sep = "\t")
  # Calculate contig size
  contigs$contig = as.character(contigs$contig)
  contigs = contigs[order(contigs$chromosome, contigs$start),] # Order table by chromosomal location
  contigs$contig_size = nchar(contigs$contig)
  colnames(contigs) = gsub(".txt.gz", "", colnames(contigs))
  return(contigs)
}

# Function for column reordering to generate BED like structure
reordering_columns = function(contigs) {
  contigs_bed = contigs %>%
    relocate(c(chromosome, start, end, strand), .before = contig)
  contigs_bed$contig_chromosome = contigs_bed$chromosome
  contigs_bed$contig_start = contigs_bed$start
  contigs_bed$contig_end = contigs_bed$end
  return(contigs_bed)
}

# Function for reading intron file
read_intron_file = function(intron_IDs) {
  # Read intron BED file as tab-separated table
  intron_IDs = read.table(intron_ID_file, sep = "\t")
  # Give column names to intron table
  colnames(intron_IDs) = c("intron_chromosome", "intron_start", "intron_end", "intron_id", "intron_score", "intron_strand")
  # Parse intron_id column to have intron transcript IDs.
  intron_IDs = intron_IDs %>%
    dplyr::select(-c(intron_score)) %>%
    separate(intron_id, remove = FALSE, into = c("transcript_id", "intron_number"), sep = "\\_intron_", extra = "merge")
  # Remove unnecessary strings from intron_number column name and increase the intron numbers by 1. 
  intron_IDs$intron_number = sub("\\_.*","",intron_IDs$intron_number) %>%
    as.numeric(intron_IDs$intron_number) + 1
  # Remove the chromosome scaffolds by eliminating rows that contain "_" in "intron_chromosome" column.
  intron_IDs = intron_IDs[!grepl("_", intron_IDs$intron_chromosome),]
  intron_IDs = intron_IDs %>%
    add_count(transcript_id, name = "number_of_introns")
  intron_IDs_bed = intron_IDs %>%
    relocate(c(intron_chromosome, intron_start, intron_end, intron_strand), .before = intron_id)
  return(intron_IDs_bed)
}

# Function for converting tables to GRanges objects
bed_to_granges = function(bed_file) { 
  granges_object = makeGRangesFromDataFrame(bed_file,
                                            keep.extra.columns = TRUE,
                                            seqnames.field = c("chromosome","chr"),
                                            start.field = c("start"),
                                            end.field = c("end"),
                                            strand.field = c("strand"))
  return(granges_object)
}

# Function for joininf count table and info table
join_tables = function(info_table, contig_table, join_column) {
  upregulated_contig_infos = dplyr::inner_join(x=info_table, 
                                               y=contig_table, 
                                               by=join_column)
  return(upregulated_contig_infos)
}


intron_summary = function(contig_infos, group_col1, control_group, test_group){
  # ➊ hangi sütunların ortalaması alınacak?
  cols_to_mean = c(control_group, test_group)
  cols_to_mean = cols_to_mean %>% purrr::flatten_chr()
  
  contig_infos %>%
    dplyr::group_by(!!rlang::ensym(group_col1)) %>%
    dplyr::reframe(
      chr                 = contig_chromosome[which.min(padj)],
      start               = contig_start   [which.min(padj)],
      end                 = contig_end     [which.min(padj)],
      intron_id, intron_number, strand, transcript_id,
      intron_gene_id, intron_gene_name, intron_biotype,
      number_of_contigs   = n_distinct(contig),
      average_contig_size = round(mean(contig_size)),
      min_contig_pvalue   = min(padj),
      mean_FC             = mean(log2FoldChange),
      
      dplyr::across(
        .cols  = dplyr::all_of(cols_to_mean),
        .fns   = mean,
        .names = "{.col}_mean"
      )
    ) %>%
    dplyr::distinct()
}

intron_summary_no_DE = function(contig_infos, group_col1, sample_name){
  cols_to_mean = as.character(unlist(sample_name))
  
  contig_infos %>%
    dplyr::group_by(!!rlang::ensym(group_col1)) %>%
    dplyr::reframe(
      chr                 = contig_chromosome[which.min(contig_start)],
      start               = min(as.numeric(contig_start)),
      end                 = max(as.numeric(contig_end)),
      intron_id, intron_number, strand, transcript_id,
      intron_gene_id, intron_gene_name, intron_biotype,
      number_of_contigs   = n_distinct(contig),
      average_contig_size = round(mean(contig_size)),
      
      dplyr::across(
        .cols  = dplyr::all_of(cols_to_mean),
        .fns   = mean,
        .names = "{.col}_mean"
      )
    ) %>%
    dplyr::distinct()
}


# Function for detection of spliced introns
find_spliced_introns = function(splice_ranges, summary_ranges){
  spliced_introns = as.data.frame(
    plyranges::join_overlap_intersect_within(splice_ranges,
                                                      bed_to_granges(ungroup(summary_ranges))),
    row.names = NULL) %>%
    dplyr::select(intron_id)
  return(spliced_introns)
}

annotate_from_gtf = function(gtf_file) {
  reader = if (grepl("\\.gz$", gtf_file)) paste("zcat", shQuote(gtf_file)) else paste("cat", shQuote(gtf_file))
  attr   = readLines(pipe(paste0(reader, " | awk -F'\\t' '$3==\"transcript\"{print $9}'")))
  pull   = function(x, key) stringr::str_match(x, paste0(key, ' "([^"]*)"'))[, 2]
  data.frame(
    ensembl_transcript_id   = sub("\\..*", "", pull(attr, "transcript_id")),
    ensembl_gene_id         = sub("\\..*", "", pull(attr, "gene_id")),
    external_gene_name      = pull(attr, "gene_name"),
    transcript_biotype      = pull(attr, "transcript_type"),
    transcript_is_canonical = as.integer(grepl('tag "Ensembl_canonical"', attr, fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}
