## This script runs fastp without any "Per read cutting by quality options" -- i.e. without trimmomatic settings.

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)

#---Functions----------------------------------------------------------------------------------------------------------------------------####
#' Function to make a data frame of fastQC paths
#' 
#' @param fastq_dir_paths Path names in which fastq files for QC are located
#' @param prefix_to_sample_name Specify if samples have been given a prefix to their original name
#' @trimmomatic_path Path name where trimmed files are stored.
#' 
#' @return Dataframe with file paths.

get_fastp_df <- function(fastq_dir_paths, prefix_to_sample_name = "", fastp_path){
  
  fastq_df <- 
    data_frame(fastq_paths_full = 
                 list.files(path = fastq_dir_paths, pattern = "fastq.gz", full.names = T),
               fastq_filename = 
                 fastq_paths_full %>% 
                 str_replace("/.*/", ""),
               sample_name = fastq_filename %>% 
                 str_replace(str_c("^", prefix_to_sample_name), "") %>% 
                 str_replace("_.*", ""), 
               fastq_paths_trimmed = fastp_path %>% str_c(., "/", fastq_filename) %>% str_replace("\\.fastq\\.gz|\\.fq\\.gz", "_trimmed.fastq.gz")) %>% 
    arrange(fastq_filename)
  
  return(fastq_df)
  
}

#' Function make directories for results
#' 
#' @param results_path File path for where output files will be stored
#' @param folder_name Folder name for output folder
#' 
#' @return Creates a new folder for output (if it doesn't already exist)

make_results_dir <- function(results_path, folder_name){
  
  results_dir_path <- str_c(results_path, "/", folder_name)
  
  if(!dir.exists(results_dir_path)){
    
    dir.create(results_dir_path)
    
  }else {
    
    print(str_c(results_dir_path, " directory already exists.."))
    
  }
  
  return(results_dir_path)
  
}

#---Main---------------------------------------------------------------------------------------------------------------------------------####

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3 && length(args) != 4) {
  
  stop("At least the first 3 arguments must be supplied:\n
       1) Directory paths to fastq separated by ','\n
       2) Output directory, where fastp and multiqc dirs will be created\n
       3) Experiment name to append onto multiqc report\n
       4) If the file paths have a prefix before the file name, this needs to be given here, otherwise defaults to empty string\n")
  
} 

# # Comment in if want to test run script
# args <- list()
# args[[1]] <- "/home/rreynolds/data/RNASeq_TestRunQC"
# args[[2]] <- "/home/rreynolds/data/RNASeq_TestRunQC/QC"
# args[[3]] <- "testrun_fastp_noperreadcutting"
# args[[4]] <- "NG-17299_"


fastq_dir_paths <- args[[1]] %>% str_split(",") %>% unlist()
output_path <- args[[2]]
experiment_name <- args[[3]]
sample_name_prefix <- ifelse(length(args) == 4, args[[4]], "")

fastp_path <- make_results_dir(results_path = output_path, folder_name = "fastp_noperreadcutting")
multiqc_path <- make_results_dir(results_path = output_path, folder_name = "multiqc")

fastq_df <- get_fastp_df(fastq_dir_paths, prefix_to_sample_name = sample_name_prefix, fastp_path)

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
    "--detect_adapter_for_pe ", # enable adapter sequence auto-detection for PE data 
    "--adapter_sequence TACACTCTTTCCCTACACGACGCTCTTCCGATCT ", # TruSeq3-PE.fa PrefixPE/1: the adapter for read1. For PE data, this is used if R1/R2 are found not overlapped.
    "--adapter_sequence_r2 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ", # TruSeq3-PE.fa PrefixPE/2: the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. 
    
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
system(command = str_c("multiqc ", fastp_path, " -o ", multiqc_path, " --ignore fastqc_*/ --ignore fastp/", " -n ", experiment_name))


