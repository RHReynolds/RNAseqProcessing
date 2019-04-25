library(tidyverse)
library(stringr)

# Set WD ----------------------------------------------------------------------------------------------

transcriptomics_diagnostics_wd <- Sys.getenv("transcriptomics_diagnostics_wd") 
setwd(transcriptomics_diagnostics_wd)

# Functions -------------------------------------------------------------------------------------------

get_bam_df <- function(bam_dir_paths){
  
  bam_df <- 
  data_frame(bam_paths_full = 
               list.files(path = bam_dir_paths, pattern = "\\.bam$", full.names = T),
             bam_bai_paths_full = 
               list.files(path = bam_dir_paths, pattern = "\\.bam\\.bai", full.names = T),
             bam_filename = 
               bam_paths_full %>% 
               str_replace("/.*/", ""),
             sample_name = bam_filename %>% 
               str_replace("_.*", ""))
    
  return(bam_df)
  
}

make_results_dir <- function(results_path, folder_name){
  
  results_dir_path <- str_c(results_path, "/", folder_name)
  
  if(!dir.exists(results_dir_path)){
    
    dir.create(results_dir_path)
    
  }else {
    
    print(str_c(results_dir_path, " directory already exists.."))
    
  }
  
  return(results_dir_path)
  
}



# Main ------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3 && length(args) != 4) {
  
  stop("At least the first 3 arguments must be supplied:\n
       1) Directory paths to BAMs separated by ','\n
       2) Output directory, where fastqc/trimmomatic and multiqc dirs will be created\n
       3) Experiment name to append onto multiqc report\n
       4) If the file paths have a prefix before the file name, this needs to be given here, otherwise defaults to empty string\n")
  
} 

# args <- list()
# args[[1]] <-  c("/data/RNAseq_diagnostics/mito_samples/postive_controls/2018-09-11", "/data/RNAseq_diagnostics/mito_samples/postive_controls/2018-09-19")
# args[[2]] <-  "/data/RNAseq_diagnostics/mito_samples/postive_controls/QC"
# args[[3]] <-  "mito_samples_positive_controls"
# args[[4]] <- "NG-17299_"

# args <- list()
# args[[1]] <-  "/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/"
# args[[2]] <-  "/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC/"
# args[[3]] <-  "/data/references/ensembl/bed/v95/ensembl_GRCh38_v95.bed"
# args[[4]] <-  "100"

bam_dir_paths <- args[[1]] %>% str_split(",") %>% unlist()
output_dir <- args[[2]]
bed_gene_model_path <- args[[3]]
read_length <- args[[4]] 

output_path <- make_results_dir(results_path = output_dir, folder_name = "RSeQC")

bam_df <- get_bam_df(bam_dir_paths)

sample_names_uniq <- bam_df$sample_name %>% unique()

print(str_c(Sys.time(), " - Performing RSeQC for samples: ", str_c(sample_names_uniq, collapse = ", ")))

system(command = str_c("/home/dzhang/.local/bin/geneBody_coverage.py",
                       " -r ", bed_gene_model_path,
                       " -i ", "/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/",
                       " -o ", output_path, "/", "geneBody_coverage_"),
       wait = F)

for(i in seq_along(sample_names_uniq)){
  
  sample_name_to_filter <- sample_names_uniq[i]
  
  print(str_c(Sys.time(), " - : ", sample_name_to_filter))

  bam_per_sample_paths <- 
    bam_df %>% 
    filter(sample_name == sample_name_to_filter) %>% 
    .[["bam_paths_full"]]
  
  bam_bai_per_sample_paths <- 
    bam_df %>% 
    filter(sample_name == sample_name_to_filter) %>% 
    .[["bam_bai_paths_full"]]
  
  if(length(bam_per_sample_paths) != 1) { stop(str_c("number bam files for ", sample_name_to_filter, " not 1 expected")) }
  
  system(command = str_c("/home/dzhang/.local/bin/clipping_profile.py -i ", bam_per_sample_paths, ' -s "PE" -o ', 
                         output_path, "/", sample_name_to_filter, "_"), 
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/inner_distance.py -i ", bam_per_sample_paths, 
                         " -r ", bed_gene_model_path,
                         " -o ", output_path, "/", sample_name_to_filter, "_"), 
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/junction_annotation.py -i ", bam_per_sample_paths, 
                         " -r ", bed_gene_model_path,
                         " -m 20 ", # minimum intron length
                         " -o ", output_path, "/", sample_name_to_filter, "_"), 
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/junction_saturation.py -i ", bam_per_sample_paths, 
                         " -r ", bed_gene_model_path,
                         " -m 20 ", 
                         " -o ", output_path, "/", sample_name_to_filter, "_"), 
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/mismatch_profile.py -i ", bam_per_sample_paths,
                         " -l ", read_length,
                         " -o ", output_path, "/", sample_name_to_filter, "_"),
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/read_distribution.py -i ", bam_per_sample_paths, 
                         " -r ", bed_gene_model_path," > ", 
                         output_path, "/", sample_name_to_filter, "_read_distribution.txt"),
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/read_duplication.py -i ", bam_per_sample_paths,
                         " -u 20000",
                         " -o ", output_path, "/", sample_name_to_filter, "_"),
         wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/read_GC.py -i ", bam_per_sample_paths,
                         " -o ", output_path, "/", sample_name_to_filter, "_"),
         wait = F)

  system(command = str_c("/home/dzhang/.local/bin/RNA_fragment_size.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path, " > ",
                         output_path, "/", sample_name_to_filter, "_RNA_fragment_size.txt"), wait = F)
  
  system(command = str_c("/home/dzhang/.local/bin/bam_stat.py -i ", bam_per_sample_paths, " > ", 
                         output_path, "/", sample_name_to_filter, "_bam_stat.txt"), 
         wait = T)
  
}



      