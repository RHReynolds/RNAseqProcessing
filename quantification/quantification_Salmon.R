# Author(s): Regina H. Reynolds

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(stringr)
library(optparse)
library(RNAseqProcessing)

# Main ------------------------------------------------------------------------------------------------

arguments <- parse_args(OptionParser(usage = "%prog",
                                     prog = "Salmon quantification for paired-end reads",
                                     description="Script for running mapping-based Salmon quantification, with selective alignment. Required inputs:\n <fastq_dir_paths>: Directory paths to fastq separated by ','.\n <salmon_index_path>: Directory path to a 'decoy-aware' salmon transcriptome index. This is generated using Salmon tools.\n <output_path>: Path to output folder.\n",
                                     option_list=list(
                                       make_option(c("-p","--sample_prefix"), default = "", help="Call flag if the fastq paths have a prefix before the sample name (as is often added by the sequencer). E.g. Fastq file name may be NM3330_PDC05_A1A2_GM-T_S4_R1_001_QC.fastq.gz. As sample name is 'PDC05_A1A2_GM-T', will need to provide the argument 'NM3330_' or 'NM...._' if want it to be generalisable to other NM tags. [default: Empty string]"),
                                       make_option(c("-s","--sample_suffix"), default = "", help="Call flag if there is text (or regex) that needs to be excluded from the tail end of the filename to get the sample name. E.g. For M3330_PDC05_A1A2_GM-T_S4_R1_001_QC.fastq.gz, would use argument '_S.*' to remove everything after _S. [default: Empty string]"),
                                       make_option(c("-l","--library_type"), default = "A", help="Call flag to supply library type. Salmon NEEDS a description of the sequencing library. As a default, Salmon infers the library type based on how the first few thousand reads map to the transcriptome. However, it is best to enter the library type based on prior knowledge. [default: A]"),
                                       make_option(c("-g", "--gene_map"), default = NULL, help = "Call flag to enable quantification at gene level. If enabled, must supply a file containing a mapping of transcripts to genes. Salmon will then output both quant.sf and quant.genes.sf files, where the latter contains aggregated gene-level abundance estimates. The transcript to gene mapping should be provided as either a GTF file, or in a simple tab-delimited format where each line contains the name of a transcript and the gene to which belongs separated by a tab. In GTF/GFF format, the 'transcript_id' is assumed to contain the transcript identifier and the 'gene_id' is assumed to contain the corresponding gene identifier. [default: NULL]")
                                       )),
                        positional_arguments = 3)


# # Comment in if want to test run script arguments
# arguments <- list()
# arguments$args[1] <- "/data/RNAseq_PD/test/QC/fastp"
# arguments$args[2] <- "/tools/salmon/salmonReferences/ensembl_v97/Homo_sapiens.GRCh38.97.cdna.all.ncrna_index"
# arguments$args[3] <- "/data/RNAseq_PD/test/salmon_quant"
# arguments$opt$sample_prefix <- "NM...._"
# arguments$opt$sample_suffix <- "_S.*"
# arguments$opt$library_type  <- "ISF"
# arguments$opt$gene_map  <- "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf"

opt <- arguments$opt

# Positional arguments
fastq_dir_paths <- arguments$args[1] %>% str_split(",") %>% unlist()
salmon_index_path <- arguments$args[2]
output_path <- arguments$args[3]

# If '-l' flag enabled, change library type accordingly.
if(!is.null(opt$library_type)){
  cat("Library type: ", opt$library_type, "\n")
  library_type <- opt$library_type
} else {
  cat("No library type provided. Salmon will infer automatically.\n")
  library_type <- "A"
}

# Create fastq dataframe
fastq_df <- RNAseqProcessing::get_fastqc_for_STAR_df(fastq_dir_paths,
                                                     prefix_to_sample_name = opt$sample_prefix,
                                                     to_exclude_to_get_sample_name = opt$sample_suffix)

sample_names_uniq <- fastq_df$sample_name %>% unique()
threads <- 32

print(str_c(Sys.time(), " - Salmon quantification for samples: ", str_c(sample_names_uniq, collapse = ", ")))

# Loop calling STAR
for(i in seq_along(sample_names_uniq)){

  sample_name_to_filter <- sample_names_uniq[i]

  print(str_c(Sys.time(), " - quantifying sample: ", sample_name_to_filter))

  fastq_per_sample_paths_trimmed_paired <-
    fastq_df %>%
    dplyr::filter(sample_name == sample_name_to_filter) %>%
    .[["fastq_paths_trimmed_paired"]]

  if(length(fastq_per_sample_paths_trimmed_paired) != 2) { stop(str_c("number fastq files for ", sample_name_to_filter, " not 2, expected because of paired-end")) }

  # Create results directory for salmon output
  results_path <- RNAseqProcessing::make_results_dir(output_path, sample_name_to_filter)

  salmon_cmd <- str_c("/tools/salmon/salmon-latest_linux_x86_64/bin/salmon quant",
                    " --libType ", library_type, # Default option (A) will make Salmon infer library type from first few thousand reads. If known, be sure to supply.
                    " --index ", salmon_index_path, # Salmon index
                    " --mates1 ", fastq_per_sample_paths_trimmed_paired[1], # Read 1
                    " --mates2 ", fastq_per_sample_paths_trimmed_paired[2], # Read 2
                    " --output ", results_path, # Salmon outputs files with default name, therefore to avoid overwriting files in loop, must save output to sample-specific directory.
                    " --threads ", threads # Number of threads to use concurrently. Salmon is designed to work well with many threads, so, if you have a sufficient number of processors, larger values here can speed up the run substantially.
                    )


  # If geneMap flag enabled
  if(!is.null(opt$gene_map)){

    # Just contains additional --geneMap argument
    salmon_cmd <- str_c(salmon_cmd ,
                      " --geneMap ", opt$gene_map  # path to GTF containing transcript mapping to genes
    )

  }

  salmon_cmd <- str_c(salmon_cmd,
                      RNAseqProcessing::get_salmon_parameters_set())

  # print(salmon_cmd)

  system(command = salmon_cmd)

}



