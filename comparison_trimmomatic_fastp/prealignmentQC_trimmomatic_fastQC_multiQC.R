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

get_fastqc_df <- function(fastq_dir_paths, prefix_to_sample_name = "", trimmomatic_path){
  
  fastq_df <- 
    data_frame(fastq_paths_full = 
                 list.files(path = fastq_dir_paths, pattern = "fastq.gz", full.names = T),
               fastq_filename = 
                 fastq_paths_full %>% 
                 str_replace("/.*/", ""),
               sample_name = fastq_filename %>% 
                 str_replace(str_c("^", prefix_to_sample_name), "") %>% 
                 str_replace("_.*", ""), 
               fastq_paths_trimmed = trimmomatic_path %>% str_c(., "/", fastq_filename), 
               fastq_paths_trimmed_paired = fastq_paths_trimmed %>% str_replace("\\.fastq\\.gz|\\.fq\\.gz", "_paired.fastq.gz"), 
               fastq_paths_trimmed_unpaired = fastq_paths_trimmed %>% str_replace("\\.fastq\\.gz|\\.fq\\.gz", "_unpaired.fastq.gz")) %>% 
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
       2) Output directory, where fastqc/trimmomatic and multiqc dirs will be created\n
       3) Experiment name to append onto multiqc report\n
       4) If the file paths have a prefix before the file name, this needs to be given here, otherwise defaults to empty string\n")
  
} 

# # Comment in if want to test run script
# args <- list()
# args[[1]] <-  "/home/rreynolds/data/RNASeq_TestRunQC"
# args[[2]] <-  "/home/rreynolds/data/RNASeq_TestRunQC/QC"
# args[[3]] <-  "testrun_fastQC_trimmomatic"
# args[[4]] <- "NG-17299_"

fastq_dir_paths <- args[[1]] %>% str_split(",") %>% unlist()
output_path <- args[[2]]
experiment_name <- args[[3]]
sample_name_prefix <- ifelse(length(args) == 4, args[[4]], "")

fastqc_pretrim_path <- make_results_dir(results_path = output_path, folder_name = "fastqc_pretrim")
trimmomatic_path <- make_results_dir(results_path = output_path, folder_name = "trimmomatic")
fastqc_posttrim_path <- make_results_dir(results_path = output_path, folder_name = "fastqc_posttrim")
multiqc_path <- make_results_dir(results_path = output_path, folder_name = "multiqc")

fastq_df <- get_fastqc_df(fastq_dir_paths, prefix_to_sample_name = sample_name_prefix, trimmomatic_path)

sample_names_uniq <- fastq_df$sample_name %>% unique()

print(str_c(Sys.time(), " - Performing QC for samples: ", str_c(sample_names_uniq, collapse = ", ")))

# run fastqc on the samples that have not yet been trimmed
# -o specifies output directory
# -t specifies number of threads to be used
system(command = str_c("fastqc ", str_c(fastq_df$fastq_paths_full, collapse = " "), " -o ", fastqc_pretrim_path, " -t 15"))

# trimming sequences using trimmomatic
for(i in seq_along(sample_names_uniq)){
  
  sample_name_to_filter <- sample_names_uniq[i]
  
  print(str_c(Sys.time(), " - Trimming: ", sample_name_to_filter))
  
  fastq_per_sample_paths <- 
    fastq_df %>% 
    filter(sample_name == sample_name_to_filter) %>% 
    .[["fastq_paths_full"]]
  
  fastq_per_sample_paths_trimmed_paired <- 
    fastq_df %>% 
    filter(sample_name == sample_name_to_filter) %>% 
    .[["fastq_paths_trimmed_paired"]] 
  
  fastq_per_sample_paths_trimmed_unpaired <- 
    fastq_df %>% 
    filter(sample_name == sample_name_to_filter) %>% 
    .[["fastq_paths_trimmed_unpaired"]] 
  
  if(length(fastq_per_sample_paths) != 2) { stop(str_c("number fastq files for ", sample_name_to_filter, " not 2, expected because of paired-end")) }
  
  # PE specifies that reads are paired end
  # need to specify the number of threads to use (depends on server)
  # need to specify the adaptor sequences that should be trimmed
  system(command = str_c("java -jar /tools/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 15 ",
                         fastq_per_sample_paths[1], " ", fastq_per_sample_paths[2], " ",
                         fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_unpaired[1], " ",  
                         fastq_per_sample_paths_trimmed_paired[2], " ", fastq_per_sample_paths_trimmed_unpaired[2], " ",
                         "ILLUMINACLIP:/tools/Trimmomatic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"))
  
}

system(command = str_c("fastqc ", str_c(str_c(fastq_df$fastq_paths_trimmed_paired, collapse = " "), " ", str_c(fastq_df$fastq_paths_trimmed_unpaired, collapse = " ")),
                       " -o ", fastqc_posttrim_path, " -t 15"))

# Assemble fastqc results in a multiqc output
system(command = str_c("multiqc ", output_path, " -o ", multiqc_path, " -n ", experiment_name))

