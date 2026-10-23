source("functions.R")

# -------------------------------- Load required libraries --------------------------------

suppressWarnings(suppressPackageStartupMessages(library(dplyr, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(stringr, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(plyranges, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(data.table, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(biomaRt, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(dbplyr, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(DESeq2, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(purrr, quietly = T)))
suppressWarnings(suppressPackageStartupMessages(library(data.table, quietly = T)))


# -------------------------------- Input and Reference files --------------------------------

# Read command line arguments
args = commandArgs(TRUE) # Enables to have command line argument

# Import input parameters
contig_file = args[1] # Input annotation table
intron_ID_file = args[2] # intron table from Nextflow parameters
metafile = args[3]

control = args[4]
test = args[5]

merged_splices = args[6]
retained_intron_file = args[7]
extended_exons_file = args[8]
gtf_file = args[9]

# -------------------------------- Convert contig annot table to BED like --------------------------------

# Conversion of the data frame into tibble to filter via dplyr
contigs <- as_tibble(read_contig_file(contig_file))

# Filter the intronic contigs and remove contigs table
intronic_contigs = dplyr::filter(contigs, is_intronic == "1")

# Separate the counts and annotation infos of contigs
col1 = "tag"
col2 = "mapped_to"
start_index = which(names(intronic_contigs) == col1)
end_index = which(names(intronic_contigs) == col2)
# contig_counts = intronic_contigs %>%
#   dplyr::select(c(names(intronic_contigs)[(start_index):(end_index - 1)], "chromosome", "start", "end", "strand"))

contig_infos = intronic_contigs %>%
  dplyr::select(-c(names(intronic_contigs)[(start_index + 1):(end_index - 1)]))

# --------------------------------------- Prepare Intron and Contig Table -------------------------------------

# Execute the read_intron_file function to have intron table
intron_IDs <- as_tibble(read_intron_file(intron_ID_file))

# ---- Gene annotation straight from the GENCODE GTF (offline; no BioMart) ----
intron_names = annotate_from_gtf(gtf_file)
intron_IDs$ids = gsub("\\..*", "", intron_IDs$transcript_id)
matched_intron_ids = match(intron_IDs$ids, intron_names$ensembl_transcript_id)
intron_IDs$intron_gene_id   = intron_names$ensembl_gene_id[matched_intron_ids]
intron_IDs$intron_gene_name = intron_names$external_gene_name[matched_intron_ids]
intron_IDs$is_canonical     = intron_names$transcript_is_canonical[matched_intron_ids]
intron_IDs$intron_size      = as.numeric(intron_IDs$intron_end) - as.numeric(intron_IDs$intron_start)
intron_IDs$intron_biotype   = intron_names$transcript_biotype[matched_intron_ids]

# Extract the canonical introns
canonical_introns = dplyr::filter(intron_IDs, is_canonical == 1)
canonical_introns_gr = bed_to_granges(canonical_introns)
intron_ids_gr = bed_to_granges(intron_IDs)

# Convert contig file to BED like table
intron_bed = reordering_columns(contig_infos)

# Convert intronic contigs table to GRanges object
intronic_contigs_gr = bed_to_granges(intron_bed)

# --------------------------------- Add Intron IDs to contigs --------------------------------------

# Give intron IDs to contigs
joined_introns = as.data.frame(plyranges::join_overlap_inner_directed(intronic_contigs_gr, canonical_introns_gr))

# Find contigs that are lost as they are aligned to noncanonical introns
lost_contigs = setdiff(intron_bed$contig, joined_introns$contig)
lost_contigs_infos = intron_bed[intron_bed$contig %in% lost_contigs,]
lost_contigs_infos_gr = bed_to_granges(lost_contigs_infos)
lost_contigs_ids = plyranges::join_overlap_inner_directed(lost_contigs_infos_gr, intron_ids_gr)
lost_contigs_ids_df = as.data.frame(lost_contigs_ids)

unique_lost_contigs = lost_contigs_ids_df %>%
  group_by(contig_chromosome, contig_start, contig_end, strand) %>%
  slice(1) %>%       # Keep only the first occurrence in each group
  ungroup() 

# Add lost contigs joined_contigs table
all_intronic_contigs = dplyr::full_join(x=joined_introns, y=unique_lost_contigs, by=colnames(joined_introns)) %>%
  dplyr::mutate(is_canonical = ifelse(intron_id %in% canonical_introns$intron_id, 1, 0)) %>%
  group_by(contig_chromosome, contig_start, contig_end, strand) %>%
  filter(
    # If there is at least one row with is_canonical == 1, keep the first one
    if (any(is_canonical == 1)) {
      row_number() == min(which(is_canonical == 1))  # Keep only the first occurrence of is_canonical == 1 within the group
    } else {
      # If no row has is_canonical == 1, keep the first row with is_canonical == 0
      row_number() == 1
    }
  ) %>%
  ungroup()

# Further filtering of the intronic contigs
# joined_introns = dplyr::filter(joined_introns, ci_ratio < 0.5 & nb_hit < 4 & contig_size > 50)
all_intronic_contigs = dplyr::filter(all_intronic_contigs, nb_hit < 4)
all_intronic_contigs = dplyr::filter(all_intronic_contigs, contig_size > 50)
all_intronic_contigs = dplyr::filter(all_intronic_contigs, is_exonic == '0')

# DESeq2 preparation

metadata = read.table(metafile, header = T, sep = "\t", check.names = F)

metadata = metadata %>% mutate(sample_name = paste(condition, replicate, sep = "_"))
row.names(metadata) = metadata$sample_name
rownames(metadata) <- sub("^X", "", rownames(metadata))


tmp =  unique(intronic_contigs %>%
  dplyr::select(c(names(intronic_contigs)[(start_index):(end_index - 1)])))

row.names(tmp) = as.character(tmp$tag)
contig_counts_for_deseq2 = tmp %>% dplyr::select(-tag)
contig_counts_for_deseq2[] <- lapply(contig_counts_for_deseq2, round, digits = 0)
row.names(contig_counts_for_deseq2) = as.character(tmp$tag)

# Reorder sampleID's to match featureCounts column order.
colnames(contig_counts_for_deseq2) <- gsub("\\.", "-", colnames(contig_counts_for_deseq2))
colnames(contig_counts_for_deseq2) <-
  sub("^X", "", colnames(contig_counts_for_deseq2))
metadata = metadata[match(colnames(contig_counts_for_deseq2), metadata$sample_name), ]

ddsMat <- DESeqDataSetFromMatrix(countData = contig_counts_for_deseq2,
                                 colData = metadata,
                                 design = ~condition)

ddsMat <- DESeq(ddsMat, quiet = T)
unique_conditions = unique(metadata$condition)

replicates = metadata %>%
  group_by(condition) %>%
  pull(sample_name)

conditions_list <- lapply(split(replicates, metadata$condition), unique)

test_index = which(names(conditions_list) == test)
control_index = which(names(conditions_list) == control)

control_group = conditions_list[control_index]
test_group = conditions_list[test_index]

stopifnot(test %in% unique_conditions, control %in% unique_conditions)

res <- results(ddsMat, 
               pAdjustMethod = "fdr", 
               alpha = 0.05, 
               contrast = c("condition", test, control))                

DE_contigs = as.data.frame(res)

if (!all(rownames(contig_counts_for_deseq2) == rownames(DE_contigs))) {
  stop("Row names of the two tables do not match!")
}

contig_counts_for_deseq2$tag = row.names(contig_counts_for_deseq2)
merged_counts_DE = cbind(contig_counts_for_deseq2, DE_contigs[, c("baseMean", "log2FoldChange", "pvalue", "padj")]) %>%
  dplyr::mutate(padj = dplyr::if_else(is.na(padj), 1, padj))

# Join count table and info tables
merged_contig_infos = join_tables(all_intronic_contigs,
                                  merged_counts_DE,
                                  join_column = "tag")

# # Mean contig counts for summary function
# mean_structure = generate_mean_structure(control_group, test_group)

# Generate summary for each group

intron_id_summary = intron_summary(merged_contig_infos, intron_id, control_group, test_group)

# Prepare splice junctions
splice_junctions = read.table(merged_splices, header = F, sep = "\t")
colnames(splice_junctions) = c("chromosome", "start", "end", "strand")
splice_junctions = splice_junctions[grepl("chr", splice_junctions$chromosome),]
splice_junctions$strand = gsub("0","*",splice_junctions$strand)
splice_junctions_ranges = bed_to_granges(splice_junctions)
colnames(splice_junctions) = c("intron_chromosome", "intron_start", "intron_end", "intron_strand")

# Prepare retained introns 
retained_introns = read.table(retained_intron_file, header = F)
colnames(retained_introns) = "intron_id"

# Prepare exon extensions
extended_exons = read.table(extended_exons_file, header = F)
colnames(extended_exons) = c("chromosome", "start", "end", "id", "score", "strand")
extended_exons_gr = bed_to_granges(extended_exons)
overlaps = findOverlaps(intron_ids_gr, extended_exons_gr)
extension_containing_introns = as.data.frame(mcols(intron_ids_gr)$intron_id[queryHits(overlaps)])
colnames(extension_containing_introns) = "intron_id"


# Find condition specific spliced introns
all_junctions = find_spliced_introns(splice_junctions_ranges, intron_IDs)

splice_containing_introns = intron_id_summary %>%
  dplyr::mutate(is_spliced = ifelse(intron_id %in% all_junctions$intron_id, 1, 0)) %>%
  dplyr::filter(is_spliced == "1")

# Remove chrM entries from splice_junctions
splice_junctions <- splice_junctions[splice_junctions$intron_chromosome != "chrM", ]

# Split splice_junctions by chromosome
splice_junctions_by_chr <- split(splice_junctions, splice_junctions$intron_chromosome)

# Split intron_IDs by chromosome
intron_IDs_by_chr <- split(intron_IDs, intron_IDs$intron_chromosome)

# Define chunk sizes
chunk_size_splice <- 1e4
chunk_size_intron <- 1e5

all_results <- rbindlist(lapply(names(splice_junctions_by_chr), function(chr) {
  # Grab chromosome-specific data
  chr_splice <- splice_junctions_by_chr[[chr]]
  chr_intron <- intron_IDs_by_chr[[chr]]
  
  # If either data frame is empty, skip
  if (nrow(chr_splice) == 0 || nrow(chr_intron) == 0) {
    return(NULL)
  }
  
  # Convert to data.table and set keys
  dt_splice <- as.data.table(chr_splice)
  dt_intron <- as.data.table(chr_intron)
  setkey(dt_splice, intron_chromosome)
  setkey(dt_intron, intron_chromosome)
  
  # Define join conditions as named character vectors
  merge_conditions <- list(
    c("intron_chromosome", "intron_start" = "intron_start"),
    c("intron_chromosome", "intron_start" = "intron_end"),
    c("intron_chromosome", "intron_end" = "intron_start"),
    c("intron_chromosome", "intron_end" = "intron_end")
  )
  
  # Perform targeted joins using .EACHI
  merged_results <- rbindlist(lapply(merge_conditions, function(cond) {
    dt_splice[dt_intron, on = cond, nomatch = 0, allow.cartesian = TRUE]
  }), use.names = TRUE, fill = TRUE)
  
  return(merged_results)
}), use.names = TRUE, fill = TRUE)


# This final data frame contains the results for all chromosomes.
possible_cassette_exons <- all_results

# Filter summary tables by retained introns and spliced introns
filtered_candidates = intron_id_summary %>%
  dplyr::mutate(
    is_spliced = ifelse(intron_id %in% all_junctions$intron_id, 1, 0),
    is_extended = ifelse(intron_id %in% extension_containing_introns$intron_id, 1, 0),
    is_retained = ifelse(intron_id %in% retained_introns$intron_id, 1, 0),
    is_casette = ifelse(intron_id %in% possible_cassette_exons$intron_id, 1, 0),
    is_canonical = ifelse(transcript_id %in% canonical_introns$transcript_id, 1,0)
  ) %>%
  dplyr::filter(is_retained == 0 & is_casette == 0 & is_extended == 0) %>%
  dplyr::rename_all(~gsub(" =.*", "", .)) %>%
  dplyr::select(-c("is_extended", "is_retained", "is_casette"))

# Write final output tables
write.table(filtered_candidates, quote = F, sep = '\t', row.names = F, col.names = T, file = paste0(test, "_vs_", control, "_candidates", ".tsv"))
save(control, test, file = "conditions.Rdata")
write.table(extension_containing_introns, quote = F, sep = '\t', row.names = F, col.names = T, file = 'extension_containing_introns.txt')
