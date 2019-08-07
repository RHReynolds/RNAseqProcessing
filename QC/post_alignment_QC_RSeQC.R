# Author(s): David Zhang, Regina H. Reynolds & Karishma D'Sa

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)
library(optparse)
library(RNAseqProcessing)

# Main ------------------------------------------------------------------------------------------------

arguments <- parse_args(OptionParser(usage = "%prog",
                                     prog = "Post-alignment QC",
                                     description="Script for sorting and indexing BAMs, and running post-alignment QC. Many downstream programs (including RSeQC) only take sorted BAM files, thus sorting and indexing by samtools performed prior to QC steps.\n Required input:\n <bam_dir_paths>: Directory paths to BAMs separated by ','\n <output_dir>: Output directory, where RSeQC directory with QC outputs will be created.\n <bed_gene_model_path>: Reference gene model in bed format.\n <read_length>: Read length.\n",
                                     option_list = list(
                                       make_option(c("-p","--sample_prefix"), default = "", help="If the bam paths have a prefix before the sample name (as is often added by the sequencer), this needs to be given here.\n [default: Empty string]"),
                                       make_option(c("-s","--sample_suffix"), default = "", help="The text (or regex) that needs to be excluded from the tail end of the filename to get the sample name.\n E.g. For 'PD115_BulkNuc-T_Aligned.sortedByCoord.out.bam', would use argument '_Aligned.sortedByCoord.out.bam' to remove this from the tail end of the file.\n [default: Empty string]")
                                     )),
                        positional_arguments = 4)


# # Comment in if want to test run script arguments
# arguments <- list()
# arguments$args[1] <- "/data/RNAseq_PD/test/STAR"
# arguments$args[2] <- "/data/RNAseq_PD/test/QC/"
# arguments$args[3] <- "/data/references/ensembl/bed/v97/ensembl_GRCh38_v97.bed"
# arguments$args[4] <- "100"
# arguments$opt$sample_prefix <- ""
# arguments$opt$sample_suffix <- "_Aligned.sortedByCoord.out.bam"

# Positional arguments
bam_dir_paths <- arguments$args[1] %>% str_split(",") %>% unlist()
output_dir <- arguments$args[2]
bed_gene_model_path <- arguments$args[3]
read_length <- arguments$args[4]

opt <- arguments$opt

bam_df <- RNAseqProcessing::get_bam_df(bam_dir_paths,
                                       prefix_to_sample_name = opt$sample_prefix,
                                       to_exclude_to_get_sample_name = opt$sample_suffix)

sample_names_uniq <- bam_df$sample_name %>% unique()

# Loop for sorting and indexing bam files-------
print(str_c(Sys.time(), " - Sorting and indexing the bams for: ", str_c(sample_names_uniq, collapse = ", ")))

for(i in seq_along(sample_names_uniq)){

  sample_name_to_filter <- sample_names_uniq[i]

  print(str_c(Sys.time(), ": ", sample_name_to_filter))

  bam_per_sample_paths <-
    bam_df %>%
    dplyr::filter(sample_name == sample_name_to_filter) %>%
    .[["bam_paths_full"]]

  if(length(bam_per_sample_paths) != 1) { stop(str_c("number bam files for ", sample_name_to_filter, " not 1 expected")) }

  ## Insert function for calling samtools sort and index -- only call one bam_dir_path

  RNAseqProcessing::call_samtools_sort_index(bam_per_sample_paths = bam_per_sample_paths,
                                             output_path = bam_dir_paths[1],
                                             sample_name = sample_name_to_filter)

  # Keep only sorted .bam file to save space and remove unsorted .bam
  print(str_c(Sys.time(), " - Removing unsorted file: ", bam_per_sample_paths))

  file.remove(bam_per_sample_paths)

}

# Loop for running RSeQC-------

output_path <- RNAseqProcessing::make_results_dir(results_path = output_dir, folder_name = "RSeQC")

# Overwriting bam_df previously used for sorting and indexing bams
# Filtering for only those bams that have been sorted and indexed
bam_df <- RNAseqProcessing::get_bam_df(bam_dir_paths,
                                       prefix_to_sample_name = opt$sample_prefix,
                                       to_exclude_to_get_sample_name = "_Aligned.sortedBysamtools.out.bam") %>%
  dplyr::filter(str_detect(bam_paths_full, pattern = "Aligned.sortedBysamtools.out.bam"))

sample_names_uniq <- bam_df$sample_name %>% unique()

print(str_c(Sys.time(), " - Performing RSeQC for samples: ", str_c(sample_names_uniq, collapse = ", ")))

# geneBody_coverage takes as input more than one .bam file, and produces a combined output for all inputted .bam files
# wait = F, as process does not require much memory + RAM, but does require a lot of time, therefore allow to run while running additional modules
system(command = str_c("/tools/RSeQC-2.6.4/scripts/geneBody_coverage.py",
                       " -r ", bed_gene_model_path,
                       " -i ", str_c(bam_df$bam_paths_full, collapse = ","),
                       " -o ", output_path, "/", "geneBody_coverage"),
       wait = F)

# Remaining RSeQC modules chosen run on individual .bam files, therefore loop through individual .bam files
# Last system command in the loop set to wait = T to prevent
for(i in seq_along(sample_names_uniq)){

  sample_name_to_filter <- sample_names_uniq[i]

  print(str_c(Sys.time(), ": ", sample_name_to_filter))

  bam_per_sample_paths <-
    bam_df %>%
    dplyr::filter(sample_name == sample_name_to_filter) %>%
    .[["bam_paths_full"]]

  if(length(bam_per_sample_paths) != 1) { stop(str_c("number bam files for ", sample_name_to_filter, " not 1 expected")) }

  # Somewhat dirty fix to ensure that process does not claim all available memory.
  # Determine available memory when running each sample
  mem_avail <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern=TRUE))
  mem_total <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern=TRUE))

  ratio_memavail_memtotal <- mem_avail/mem_total * 100

  if(ratio_memavail_memtotal <= 50){

    # A few more of the system commands have wait = T, which means process must finish before next system command called
    print(str_c("Memory available: ", round(ratio_memavail_memtotal), "%. Therefore, running memory conserving script."))

    RNAseqProcessing::call_RSeQC_modules(bam_per_sample_paths = bam_per_sample_paths,
                                         bed_gene_model_path = bed_gene_model_path,
                                         output_path = output_path,
                                         read_length = read_length,
                                         sample_name = sample_name_to_filter,
                                         wait_flag = TRUE)

  } else{

    # Only wait = T is for the last process
    print(str_c("Memory available: ", round(ratio_memavail_memtotal), "%. Therefore, running script with normal memory options."))

    RNAseqProcessing::call_RSeQC_modules(bam_per_sample_paths = bam_per_sample_paths,
                                         bed_gene_model_path = bed_gene_model_path,
                                         output_path = output_path,
                                         read_length = read_length,
                                         sample_name = sample_name_to_filter)

  }

}


