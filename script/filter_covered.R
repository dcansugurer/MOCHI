suppressWarnings(suppressPackageStartupMessages(library(dplyr, quietly = T)))

args = commandArgs(TRUE)
load("conditions.Rdata")

# Import input parameters
filtered_intron_file = args[1] # Input annotation table
not_covered_intron_file = args[2]

filtered_intron = read.table(filtered_intron_file, sep = "\t", header = T)
not_covered_intron = read.table(not_covered_intron_file, sep = "\t", header = F)
colnames(not_covered_intron) = "intron_id"

# Filter summary tables by retained introns and spliced introns
covered_intron = filtered_intron %>%
  dplyr::mutate(
    is_covered = ifelse(intron_id %in% not_covered_intron$intron_id, 0, 1)
  ) %>%
  dplyr::filter(is_covered == 1) %>%
  dplyr::select(-c("is_covered"))


out_label = if (exists("sample_name")) sample_name else paste0(test, "_vs_", control)
write.table(covered_intron, quote = F, sep = '\t', row.names = F, col.names = T,
            file = paste0(out_label, "_covered_candidates", ".tsv"))
#write.table(covered_intron, quote = F, sep = '\t', row.names = F, col.names = T, file = paste0(test, "_vs_", control, "_covered_candidates", ".tsv"))

