# Author(s): Regina H. Reynolds

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(optparse)
library(tidyverse)

# Main ------------------------------------------------------------------------------------------------

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file_dir",
                                     description="\n Wrap around for LeafCutter differential splicing command line tool that allows multiple pairwise comparisons. \n \n Required inputs:\n <counts_file>: Intron usage counts file. Must be .txt or .txt.gz, output from clustering pipeline.\n <groups_file_dir>: Path to directory containing LeafCutter appropriate _group_files.txt with various pairwise combinations to be performed with the same count file. For each pairwise comparison, a group file must be present. Each group file should be end with the suffix '_group_file.txt' and should be a two+K column file: \n 1. sample names (must match column names in counts_file), \n 2. groups (currently only two groups, i.e. pairwise, supported. \n Some samples in counts_file can be missing from this file, in which case they will not be included in the analysis. Additional columns can be used to specify confounders, e.g. batch/sex/age. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you a categorical variable.",
                                     option_list=list(
                                       make_option(c("-o","--output_prefix"), default = "leafcutter_ds", help="The prefix for the two output files, <prefix>_cluster_significance.txt (containing test status, log likelihood ratio, degree of freedom, and p-value for each cluster) and <prefix>_effect_sizes.txt (containing the effect sizes for each intron)  [default %default]"),
                                       make_option(c("-s","--max_cluster_size"), default=Inf, help="Don't test clusters with more introns than this [default %default]"),
                                       make_option(c("-i","--min_samples_per_intron"), default=5, help="Ignore introns used (i.e. at least one supporting read) in fewer than n samples [default %default]") ,
                                       make_option(c("-g","--min_samples_per_group"), default=3, help="Require this many samples in each group to have at least min_coverage reads [default %default]"),
                                       make_option(c("-c","--min_coverage"), default=20, help="Require min_samples_per_group samples in each group to have at least this many reads [default %default]"),
                                       make_option(c("-t","--timeout"), default=30, help="Maximum time (in seconds) allowed for a single optimization run [default %default]"),
                                       make_option(c("-p","--num_threads"), default=1, help="Number of threads to use [default %default]"),
                                       make_option(c("-e","--exon_file"), default=NULL, help="File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name. Optional, only just to label the clusters."))),
                        positional_arguments = 2)


# # Comment in if want to test run script arguments
# arguments <- list()
# arguments$args[1] <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/intron_clustering/tissue_polyA_test_diseasegroups_perind_numers.counts.gz"
# arguments$args[2] <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/group_files/"
# arguments$opt$output_prefix <- "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/"
# arguments$opt$max_cluster_size <- "Inf"
# arguments$opt$min_samples_per_intron  <- "5"
# arguments$opt$min_samples_per_group  <- "3"
# arguments$opt$min_coverage <- "20"
# arguments$opt$timeout <- "30"
# arguments$opt$num_threads <- "15"
# arguments$opt$exon_file <-"/data/references/ensembl/gtf_gff3/v97/leafcutter/Homo_sapiens.GRCh38.97_LC_exon_file.txt.gz"

opt <- arguments$opt

# Positional arguments
count_file <- arguments$args[1]
groups_file_dir <- arguments$args[2]

# Dataframe of _group_file.txt files in directory
group_file_df <- tibble(group_file_path = list.files(path = groups_file_dir, pattern = "_group_file.txt", full.names = T),
                        group_file_name = list.files(path = groups_file_dir, pattern = "_group_file.txt", full.names = T) %>%
                          str_replace("/.*/", "") %>%
                          str_replace("_group_file.txt", ""))

# Loop calling leafcutter
for(i in 1:nrow(group_file_df)){

  print(str_c("Performing leafcutter differential splicing using the group file entitled: ", group_file_df$group_file_name[i]))

  leafcutter_cmd <- str_c("Rscript /tools/leafcutter/scripts/leafcutter_ds.R ",
                          count_file, " ", # path to count file
                          group_file_df$group_file_path[i], # path to group file
                          " --output_prefix ", str_c(opt$output_prefix, "/", group_file_df$group_file_name[i]), # comparison-specific output prefix
                          " --max_cluster_size ", opt$max_cluster_size,
                          " --min_samples_per_intron ", opt$min_samples_per_intron,
                          " --min_samples_per_group ", opt$min_samples_per_group,
                          " --min_coverage ", opt$min_coverage,
                          " --timeout ", opt$timeout,
                          " --num_threads ", opt$num_threads,
                          " --exon_file ", opt$exon_file
                          )

  # print(leafcutter_cmd)

  system(command = leafcutter_cmd)

}



