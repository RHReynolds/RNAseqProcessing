## This script runs fastp without any "Per read cutting by quality options" -- i.e. without trimmomatic settings.
## No adapter trimming

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)
library(RNAseqProcessing)

#---Main---------------------------------------------------------------------------------------------------------------------------------####

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3 | length(args) > 5) {
  stop("At least the first 3 arguments must be supplied:\n
       1) Directory paths to fastq separated by ','\n
       2) Output directory, where fastp and multiqc dirs will be created\n
       3) Experiment name to append onto multiqc report\n
       4) If the file paths have a prefix before the file name, this needs to be given here, otherwise defaults to empty string\n
       5) The text (or regex) that needs to be excluded from the tail end of the filename to get the sample name. If no argument provided, defaults to empty string")

}

# # Comment in if want to test run script
# args <- list()
# args[[1]] <- "/data/RNAseq_PD/tissue_polyA_samples/raw_data"
# args[[2]] <- "/data/RNAseq_PD/tissue_polyA_samples/QC"
# args[[3]] <- "PD_tissue_polyA"
# args[[4]] <- "NM...._"


fastq_dir_paths <- args[[1]] %>% str_split(",") %>% unlist()
output_path <- args[[2]]
experiment_name <- args[[3]]
sample_name_prefix <- ifelse((length(args) >=4 && args[[4]] != "NA" ), args[[4]], "")
exclude_to_get_sample_name <- ifelse((length(args) == 5 && args[[5]] != "NA" ), args[[5]], "")

fastp_path <- RNAseqProcessing::make_results_dir(results_path = output_path, folder_name = "fastp")
multiqc_path <- RNAseqProcessing::make_results_dir(results_path = output_path, folder_name = "multiqc")

fastq_df <- RNAseqProcessing::get_fastp_df(fastq_dir_paths, prefix_to_sample_name = sample_name_prefix, fastp_path, to_exclude_to_get_sample_name = exclude_to_get_sample_name)

sample_names_uniq <- fastq_df$sample_name %>% unique()

print(str_c(Sys.time(), " - Performing QC and trimming for samples: ", str_c(sample_names_uniq, collapse = ", ")))

# run fastp
for(i in seq_along(sample_names_uniq)){

  sample_name_to_filter <- sample_names_uniq[i]

  print(str_c(Sys.time(), " - QC and Trimming: ", sample_name_to_filter))

  fastq_per_sample_paths <-
    fastq_df %>%
    filter(sample_name == sample_name_to_filter) %>%
    .[["fastq_paths_full"]]

  fastq_per_sample_paths_trimmed <-
    fastq_df %>%
    filter(sample_name == sample_name_to_filter) %>%
    .[["fastq_paths_trimmed"]]

  if(length(fastq_per_sample_paths) != 2) { stop(str_c("number fastq files for ", sample_name_to_filter, " not 2, expected because of paired-end")) }

  system(command = str_c(

    ## Input/output options
    "fastp --in1 ",
    fastq_per_sample_paths[1], " ", # read1 input file name
    "--out1 ",
    fastq_per_sample_paths_trimmed[1], " ", # read1 output file name
    "--in2 ",
    fastq_per_sample_paths[2], " ", # read2 input file name
    "--out2 ",
    fastq_per_sample_paths_trimmed[2], " ", # read2 output file name

    ## Adapter trimming options
    "--disable_adapter_trimming ", # disable adapter trimming

    ## Per read cutting by quality options
    ## N.B. Performing per read cutting by quality score interferes with deduplication by fastp, as algorithms rely on exact matching between coordinate regions of paired reads
    # "--cut_front ", # equivalent to Trimmomatic's LEADING option, if window size set to 1
    # "--cut_front_window_size 1 ",
    # "--cut_front_mean_quality 3 ", # Set to default trimmomatic quality for LEADING
    # "--cut_right ", # similar to Trimmomatic SLIDINGWINDOW option
    # "--cut_right_window_size 4 ",
    # "--cut_right_mean_quality 15 ", # the mean quality requirement option for cut_right. Range: 1~36 default: 20 (Q20), but set to 15 to match Trimmomatic.

    ## Quality filtering options -- enabled by default
    "--qualified_quality_phred 15 ", # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
    "--unqualified_percent_limit 40 ", # how many percents of bases are allowed to be unqualified (0~100). Default 40.
    "--n_base_limit 5 ", # if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5

    ## Length filtering options
    "--length_required 36 ", # reads shorter than length_required will be discared. Default is 15, but set to 36 to match Trimmomatic.

    ## Base correction by overlap analysis options
    "--correction ", # enables base correction
    "--overlap_len_require 30 ", # the minimum length of the overlapped region for overlap analysis based adapter trimming and correction. 30 by default.
    "--overlap_diff_limit 5 ", # the maximum difference of the overlapped region for overlap analysis based adapter trimming and correction. 5 by default.

    ## Overrepresented sequence analysis
    "--overrepresentation_analysis ", # enable overrepresented sequence analysis
    "--overrepresentation_sampling 20 ", #  one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20.

    ## Reporting options
    "--html ", fastp_path, "/", sample_name_to_filter, "_fastp.html ",
    "--json ", fastp_path, "/", sample_name_to_filter, "_fastp.json ",
    "--report_title '", sample_name_to_filter, "' ",

    ## Threading options
    "--thread 16"  # fastp uses up to 16 threads

  ))

}

# fastqc of fastp trimming
system(command = str_c("fastqc ", str_c(str_c(fastq_df$fastq_paths_trimmed, collapse = " ")),
                       " -o ", fastp_path, " -t 20"))

# multiqc to collect fastp and fastqc outputs
system(command = str_c("multiqc ", fastp_path, " -o ", multiqc_path, " -n ", experiment_name))


