# Author(s): Regina H. Reynolds

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(optparse)
library(tidyverse)

# Main ------------------------------------------------------------------------------------------------

arguments <- parse_args(OptionParser(usage="%prog [options] <name>_perind_numers.counts.gz <path>_cluster_significance.txt <path>_effect_sizes.txt annotation_code",
                                     description = "\n Wrap around function for the ´prepare_results.R´ Leafviz script that prepares differential splicing results for use with the Leafviz shiny platform. This script allows processing of multiple diffential splicing analyses, as opposed to just one. \n \n Required inputs include: \n <name>_perind_numers.counts.gz: The per individual intron count file outputted by LeafCutter intron clustering pipeline. \n <path>_cluster_significance.txt: Path to the directory containing the _cluster_significance.txt files outputted by LeafCutter differential splicing analyses. \n <path>_effect_sizes.txt: Path to the directory containing the _effect_sizes.txt files outputted by LeafCutter differential splicing analyses. \n annotation_code: The prefix for the annotation files that LeafCutter should use to annotate junctions to genes. The annotation_code should be something like annotation_codes/gencode_hg19/gencode_hg19.",
                                     option_list=list(
                                       make_option(c("-o", "--output_dir"), default=NULL, help = "The output directory for leafviz prepare_results."),
                                       make_option(c("--group_file_dir"), default=NULL, help="The support files used in the differential splicing analysis. Columns should be file name and condition"),
                                       make_option(c("-f","--FDR"), default=0.05, help = "the adjusted p value threshold to use [default %default]"))),
                        positional_arguments = 4)


# # Comment in if want to test run script arguments
# arguments <- list()
# arguments$args[1] <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/intron_clustering/tissue_polyA_test_diseasegroups_perind_numers.counts.gz"
# arguments$args[2] <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/"
# arguments$args[3] <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/"
# arguments$args[4] <- "/data/references/ensembl/gtf_gff3/v97/leafcutter/"
# arguments$opt$output_dir <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/leafviz/"
# arguments$opt$group_file_dir <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/group_files/"
# arguments$opt$FDR  <- "0.05"

opt <- arguments$opt

# Positional arguments
count_file <- arguments$args[1]
cluster_significant_dir <- arguments$args[2]
effect_size_dir <- arguments$args[3]
annotation_code <- arguments$args[4]

# Dataframe of ds runs
run_df <- tibble(cluster_significant_file = list.files(path = cluster_significant_dir, pattern = "_cluster_significance.txt", full.names = T),
                 effect_size_file = list.files(path = effect_size_dir, pattern = "_effect_sizes.txt", full.names = T),
                 group_file = list.files(opt$group_file_dir, pattern = "_group_file.txt", full.names = T),
                 comparison_name = list.files(path = cluster_significant_dir, pattern = "_cluster_significance.txt", full.names = T) %>%
                   str_replace("/.*/", "") %>%
                   str_replace("_cluster_significance.txt", ""))

# Loop calling prepare_results for leafviz
for(i in 1:nrow(run_df)){

  print(str_c("Performing leafviz prepare_results for: ", run_df$comparison_name[i]))

  leafviz_cmd <- str_c("Rscript /tools/leafcutter/leafviz/prepare_results.R ",
                       count_file, " ",
                       run_df$cluster_significant_file[i], " ",
                       run_df$effect_size_file[i], " ",
                       annotation_code,
                       " --output ", opt$output_dir, "/", run_df$comparison_name[i], ".Rda ", # comparison-specific output prefix
                       " --meta_data_file ", run_df$group_file[i],
                       " --FDR ", opt$FDR,
                       " --code ", run_df$comparison_name[i]
  )

  # print(leafviz_cmd)

  system(command = leafviz_cmd)

}



