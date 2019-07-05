# Author: Regina H. Reynolds & David Zhang

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)
library(RNAseqProcessing)

# Main ------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3 && length(args) != 4) {

  stop("At least the first 3 arguments must be supplied:\n
       1) Directory paths to fastq separated by ','\n
       2) Path to genome index used by STAR for mapping reads\n
       3) Path to output folder\n
       4) If the file paths have a prefix before the file name, this needs to be given here, otherwise defaults to empty string\n
       5) The text (or regex) that needs to be excluded from the tail end of the filename to get the sample name. If no argument provided, defaults to empty string\n")

}

# # Comment in if want to test run script
# args <- list()
# args[[1]] <- "/data/RNAseq_PD/tissue_polyA_samples/QC/fastp"
# args[[2]] <- "/data/STAR_data/genome_index_hg38"
# args[[3]] <- "/data/RNAseq_PD/tissue_polyA_samples/STAR"
# args[[4]] <- "NM...._"
# args[[5]] <- "_.*"

# args <- list()
# args[[1]] <- "/data/RNAseq_PD/test/QC/fastp"
# args[[2]] <- "/data/STAR_data/genome_index_hg38"
# args[[3]] <- "/data/RNAseq_PD/test/STAR"
# args[[4]] <- "NM...._"
# args[[5]] <- "_.*"

fastq_dir_paths <- args[[1]] %>% str_split(",") %>% unlist()
genome_index_path <- args[[2]]
output_path <- args[[3]]
sample_name_prefix <- ifelse((length(args) >=4 && args[[4]] != "NA" ), args[[4]], "")
exclude_to_get_sample_name <- ifelse((length(args) == 5 && args[[5]] != "NA" ), args[[5]], "")

fastq_df <- RNAseqProcessing::get_fastqc_for_STAR_df(fastq_dir_paths,
                                                     prefix_to_sample_name = sample_name_prefix,
                                                     to_exclude_to_get_sample_name = exclude_to_get_sample_name)

sample_names_uniq <- fastq_df$sample_name %>% unique()
# threads_STAR <- floor(50/length(sample_names_uniq))
# threads_STAR <- ifelse(threads_STAR > 15, 15, threads_STAR) # allow for maximum number of 15 threads per sample
threads_STAR <- 15

print(str_c(Sys.time(), " - STAR alignment for samples: ", str_c(sample_names_uniq, collapse = ", ")))

for(i in seq_along(sample_names_uniq)){

  sample_name_to_filter <- sample_names_uniq[i]

  print(str_c(Sys.time(), " - aligning sample: ", sample_name_to_filter))

  fastq_per_sample_paths_trimmed_paired <-
    fastq_df %>%
    filter(sample_name == sample_name_to_filter) %>%
    .[["fastq_paths_trimmed_paired"]]


  # This needs modifying for use with a generalised df with the necessary RG information
  RG_ID <- fastq_df %>%
    filter(sample_name == sample_name_to_filter) %>%
    .[["fastq_filename"]] %>%
    str_replace("_R.*", "") %>%
    unique()

  if(length(fastq_per_sample_paths_trimmed_paired) != 2) { stop(str_c("number fastq files for ", sample_name_to_filter, " not 2, expected because of paired-end")) }

  system(command = str_c("STAR",
                         " --runThreadN ", threads_STAR,
                         " --genomeDir ", genome_index_path,
                         " --readFilesIn ", fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_paired[2],
                         " --readFilesCommand  zcat ", # because fastq's are zipped
                         "--outFileNamePrefix ",  str_c(output_path, "/", sample_name_to_filter, "_ "),
                         "--outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
                         "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM
                         "--outSAMattrRGline ", str_c("ID:", RG_ID, " PU:xxx ", "SM:", sample_name_to_filter, " PL:Illumina LB:xxx "), #
                         "--outFilterType BySJout ", # removes spurious split reads
                         "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
                         "--outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
                         "--outFilterMismatchNoverReadLmax 0.04 ", # max number of mismatches per pair relative to read length. As per current ENCODE options.
                         "--alignIntronMin 20 ", # min intron length. As per ENCODE options.
                         "--alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
                         "--alignMatesGapMax 1000000", # max gap between pair mates. As per ENCODE options.
                         "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
                         "--alignSJDBoverhangMin 3 " # minimum annotated split read anchor. Default is 3.
                         ))
}

