# Author(s): Regina H. Reynolds

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)
library(optparse)
library(RNAseqProcessing)

# Main ------------------------------------------------------------------------------------------------

arguments <- parse_args(OptionParser(usage = "%prog",
                                     prog = "STAR: merging SJ.out.tab",
                                     description="Script for merging SJ.out.tab files outputted by alignment, removing duplicated junctions, and outputting a single file with non-duplicated junctions from all samples.\n Required inputs:\n <sj_dir_path>: Path to directory containing splice junctions files (SJ.out.tab).\n",
                                     option_list=list(
                                       make_option(c("-o","--output_path"), default = NULL, help="If different output path desired, specify here. Default is to output the merged file to the same folder wherein SJ.out.tab files are located."))),
                        positional_arguments = 1)

# # Comment in if want to test run script
# arguments <- list()
# arguments$args[1] <- "/data/RNAseq_PD/tissue_polyA_samples/STAR"
# arguments$opt$output_path  <- NULL

# Positional arguments
sj_dir_path <- arguments$args[1]

# Optional arguments
opt <- arguments$opt

# Load SJ.tab.out files
sj_df <- RNAseqProcessing::load_sj_df(sj_dir_path = sj_dir_path)

unique_junc <- sj_df %>%
  tidyr::unite(col = "junction", chr, intron_start, intron_end, sep = ":", remove = FALSE) %>%
  dplyr::distinct(junction, .keep_all = T) %>%
  dplyr::select(chr, intron_start, intron_end)

# If '-o' flag enabled, output files to different output path
if(!is.null(opt$output_path)){
  cat("Writing file to:", opt$output_path, "\n")
  write_delim(x = unique_junc, path = str_c(output_path, "/all_samples_non_duplicated_junctions.SJ.out.tab"), col_names = F, delim = "\t")
} else {
  cat("Writing file to:", sj_dir_path, "\n")
  write_delim(x = unique_junc, path = str_c(sj_dir_path, "/all_samples_non_duplicated_junctions.SJ.out.tab"), col_names = F, delim = "\t")
}
