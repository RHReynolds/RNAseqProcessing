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
                                       make_option(c("-o","--output_path"), default = "", help="If different output path desired, specify here. Default is to output the merged file to the same folder wherein SJ.out.tab files are located."),
                                       make_option(c("-f,", "--filter_number_of_samples"), default = NULL, help = "Add filter for the number of samples a junction must be found in to be included in the final merged file. E.g. if set to 5, then only junctions found in 5 or more samples included."))),
                        positional_arguments = 1)

# # Comment in if want to test run script arguments
# arguments <- list()
# arguments$args[1] <- "/data/RNAseq_PD/tissue_polyA_samples/STAR/SJ_out_1pass"
# arguments$opt$output_path  <- ""
# arguments$opt$filter_number_of_samples <- 5

# Positional arguments
sj_dir_path <- arguments$args[1]

# Optional arguments
opt <- arguments$opt

# Load SJ.tab.out files
sj_df <- RNAseqProcessing::load_sj_df(sj_dir_path = sj_dir_path)

# If '-f' flag enabled, output files to different output path
if(!is.null(opt$filter_number_of_samples)){

  cat("Filtering for junctions appearing in >=", opt$filter_number_of_samples, "samples. \n")

  filtered_junc <- sj_df %>%
    tidyr::unite(col = "junction", chr, intron_start, intron_end, strand, sep = ":", remove = FALSE) %>%
    dplyr::group_by(junction) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n >= opt$filter_number_of_samples) %>% # filter for number of samples junction is present in
    tidyr::separate(junction, c("chr", "intron_start", "intron_end", "strand"), sep = ":") %>%
    dplyr::mutate(strand = ifelse(strand == 0, ".",
                                  ifelse(strand == 1, "+", "-"))) %>%
    dplyr::select(-n)


} else {

  cat("No filter set; no filtering occurring.\n")

  filtered_junc <- sj_df %>%
    tidyr::unite(col = "junction", chr, intron_start, intron_end, strand, sep = ":", remove = FALSE) %>%
    dplyr::distinct(junction, .keep_all = T) %>%
    dplyr::select(chr, intron_start, intron_end, strand) %>%
    dplyr::mutate(strand = ifelse(strand == 0, ".",
                                  ifelse(strand == 1, "+", "-")))

}

# If '-o' flag enabled, output files to different output path
if(!opt$output_path == ""){

  cat("Writing file to:", opt$output_path, "\n")
  write_delim(x = filtered_junc, path = str_c(opt$output_path, "/merged_junctions.SJ.out.tab"), col_names = F, delim = "\t")

} else {

  cat("Writing file to:", sj_dir_path, "\n")
  write_delim(x = filtered_junc, path = str_c(sj_dir_path, "/merged_junctions.SJ.out.tab"), col_names = F, delim = "\t")

}
