#' Function to make a data frame of fastQC paths
#'
#' @param fastq_dir_paths Path names in which fastq files for QC are located
#' @param prefix_to_sample_name Specify if samples have been given a prefix to
#'   their original name
#' @trimmomatic_path Path name where trimmed files are stored.
#' @to_exclude_to_get_sample_name Exclude this text or regex from the filename
#'   to get the sample name. To substitute for an empty string, use NA.
#'   default="_.*"
#'
#' @return Dataframe with file paths.

get_fastp_df <- function(fastq_dir_paths, prefix_to_sample_name = "", fastp_path, to_exclude_to_get_sample_name="_.*"){
  
  fastq_df <- 
    data_frame(fastq_paths_full = 
                 list.files(path = fastq_dir_paths, pattern = "fastq.gz", full.names = T),
               fastq_filename = 
                 fastq_paths_full %>% 
                 str_replace("/.*/", ""),
               sample_name = fastq_filename %>% 
                 str_replace(str_c("^", prefix_to_sample_name), "") %>% 
                 str_replace(to_exclude_to_get_sample_name, ""),
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