#' Function to return a string containing the STAR parameters --sjdbFileChrStartEnd and --limitSjdbInsertNsj for sjdb.
#'
#' @param sj_file - the file containing the splice jucntions fomr the 1st pass
#'
#' @return string with the sjdb STAR parameters
#' @export
#'

get_star_parameters_sjdb<- function(sj_file){
  # Determine limitSjdb for the --limitSjdbInsertNsj flag
  limitSjdb <- get_limitSjdb(sj_file)

  sjdbStarParameter <- str_c(" --sjdbFileChrStartEnd ", sj_file,
                             " --limitSjdbInsertNsj ", format(limitSjdb, scientific=F) # maximum number of junction to be inserted to the genome on the ﬂy at the mapping stage, including those from annotations. Default is 1,000,000 -- but may need to be larger depending on annotation file.
  )
  return(sjdbStarParameter)
}


#' Function to make a data frame of trimmed paths for STAR.
#'
#' @param fastq_dir_paths Path names in which trimmed fastq files for alignment
#'   are located.
#' @param prefix_to_sample_name Specify if samples have been given a prefix to
#'   their original name.
#' @param to_exclude_to_get_sample_name The text (or regex) that needs to be
#'   excluded from the tail end of the filename to get the sample name. If no
#'   argument provided, defaults to empty string.
#'
#' @return Dataframe with file paths.
#' @export
#'

get_fastqc_for_STAR_df <- function(fastq_dir_paths, prefix_to_sample_name = "", to_exclude_to_get_sample_name = ""){

  fastq_df <-
    data_frame(fastq_paths_trimmed_paired =
                 list.files(path = fastq_dir_paths, pattern = ".fastq.gz", full.names = T),
               fastq_filename =
                 fastq_paths_trimmed_paired %>%
                 str_replace("/.*/", ""),
               sample_name = fastq_filename %>%
                 str_replace(str_c("^", prefix_to_sample_name), "") %>%
                 str_replace(str_c(to_exclude_to_get_sample_name,"$"), "")) %>%
    arrange(fastq_filename)

  return(fastq_df)

}

#' Function to make a data frame of bam paths for sorting and indexing.
#'
#' @param bam_dir_paths Path names in which bam files are located.
#' @param prefix_to_sample_name Specify if samples have been given a prefix to
#'   their original name.
#' @param to_exclude_to_get_sample_name The text (or regex) that needs to be
#'   excluded from the tail end of the filename to get the sample name. If no
#'   argument provided, defaults to empty string.
#'
#' @return Dataframe with file paths.
#' @export

get_bam_df <- function(bam_dir_paths, prefix_to_sample_name = "", to_exclude_to_get_sample_name = ""){

  bam_df <-
    data_frame(bam_paths_full =
                 list.files(path = bam_dir_paths, pattern = "\\.bam$", full.names = T),
               bam_filename =
                 bam_paths_full %>%
                 str_replace("/.*/", ""),
               sample_name = bam_filename %>%
                 str_replace(str_c("^", prefix_to_sample_name), "") %>%
                 str_replace(str_c(to_exclude_to_get_sample_name,"$"), ""))

  return(bam_df)

}

#' Load splice junctions from first pass STAR alignment.
#'
#' @param sj_dir_path Path to directory containing splice junctions files
#'   (SJ.out.tab).
#'
#' @return Returns a dataframe of all splice junctions from SJ.out.tab files.
#' @export
#'

load_sj_df <- function(sj_dir_path){

  paths <- list.files(path = sj_dir_path, pattern = "SJ.out.tab", full.names = TRUE)

  for(i in 1:length(paths)){

    sample_name <- paths[i] %>%
      str_replace("/.*/", "") %>%
      str_replace("_SJ.out.tab", "")

    cat("Loading splice junctions from:", paths[i],"\n")
    sj_df <- read_delim(file = paths[i],
                        delim = "\t",
                        col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "in_annotation", "unique_reads_junction", "multi_map_reads_junction", "max_splice_alignment_overhang"),
                        col_types = "cdddddddd") %>%
      dplyr::mutate(Sample = sample_name)

    if(i == 1){
      master_df <- sj_df
    } else {
      master_df <- master_df %>%
        bind_rows(sj_df)
    }

  }

  return(master_df)

}

#' Get fixed STAR parameters.
#'
#' Function to return a string containing the non-sampleID-dependent STAR
#' parameters. To be used in the 1st and 2nd pass alignment.
#'
#' @return String with the fixed STAR parameters.
#' @export
#'

get_star_parameters_set <-function(){

  return(str_c(" --outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
               "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM,
               "--outFilterType BySJout ", # removes spurious split reads
               "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
               "--outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
               "--outFilterMismatchNoverReadLmax 0.04 ", # max number of mismatches per pair relative to read length. As per current ENCODE options.
               "--alignIntronMin 20 ", # min intron length. As per ENCODE options.
               "--alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
               "--alignMatesGapMax 1000000 ", # max gap between pair mates. As per ENCODE options.
               "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
               "--alignSJDBoverhangMin 3 " # minimum annotated split read anchor. Default is 3.
  ))

}

#' Call samtools sort and index tools.
#'
#' @param bam_per_sample_paths File path for an individual sample .bam file.
#' @param output_path Output path for sorted and indexed .bam file
#' @param sample_name Sample name.
#'
#' @return Runs system commands for sorting and indexing .bam file.
#' @export
#'

call_samtools_sort_index <- function(bam_per_sample_paths, output_path, sample_name){

  system(command = str_c("samtools sort -m 1000000000 ", bam_per_sample_paths, " -o ",
                         output_path, "/", sample_name, "_Aligned.sortedBysamtools.out.bam"))

  system(command = str_c("samtools index ", output_path, "/", sample_name, "_Aligned.sortedBysamtools.out.bam"))

}



#' Call RSeQC modules.
#'
#' @param bam_per_sample_paths File path for an individual sample .bam file.
#' @param bed_gene_model_path Reference gene model in bed format.
#' @param output_path Output directory, where RSeQC outputs will be outputted.
#' @param read_length Read length.
#' @param sample_name Sample name.
#' @param wait_flag Call flag with "TRUE", if memory available is <= 50%.
#'
#' @return Runs system commands for RSeQC modules.
#' @export

call_RSeQC_modules <- function(bam_per_sample_paths, bed_gene_model_path, output_path, read_length, sample_name, wait_flag = NULL){

  # Construct vector of T/F for modules, depending on memory availability.
  if(!is.null(wait_flag)){

    # If wait_flag = T, this will mean a few more of the system commands have wait = T.
    # wait = T in system() will force the system to wait for the to finish before next system command called.
    wait_flags <- c(F, F, T, F, F, T, F, F, F, T)

  } else{

    # All, but the last system command, are set to wait = F.
    # Thus, all system commands will be called simultaneously and run in parallel.
    # Only last process has wait = T, to force system to wait for this process to finish.
    wait_flags <- c(F, F, F, F, F, F, F, F, F, T)

  }

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/clipping_profile.py -i ", bam_per_sample_paths, ' -s "PE" -o ',
                         output_path, "/", sample_name),
         wait = wait_flags[1])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/inner_distance.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path,
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[2])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/junction_annotation.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path,
                         " -m 20", # minimum intron length, as per STAR settings
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[3])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/junction_saturation.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path,
                         " -m 20", # minimum intron length, as per STAR settings
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[4])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/mismatch_profile.py -i ", bam_per_sample_paths,
                         " -l ", read_length,
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[5])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/read_distribution.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path," > ",
                         output_path, "/", sample_name, "_read_distribution.txt"),
         wait = wait_flags[6])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/read_duplication.py -i ", bam_per_sample_paths,
                         " -u 20000", # upper limit of reads' occurrence. Limit only used for plotting.
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[7])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/read_GC.py -i ", bam_per_sample_paths,
                         " -o ", output_path, "/", sample_name),
         wait = wait_flags[8])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/RNA_fragment_size.py -i ", bam_per_sample_paths,
                         " -r ", bed_gene_model_path, " > ",
                         output_path, "/", sample_name, "_RNA_fragment_size.txt"),
         wait = wait_flags[9])

  system(command = str_c("/tools/RSeQC-2.6.4/scripts/bam_stat.py -i ", bam_per_sample_paths, " > ",
                         output_path, "/", sample_name, "_bam_stat.txt"),
         wait = wait_flags[10])

}


#' Get fixed Salmon parameters.
#'
#' Function to return a string containing the non-sampleID-dependent salmon
#' parameters. To be used with salmon quantification.
#'
#' @return String with the fixed salmon parameters.
#' @export
#'

get_salmon_parameters_set <-function(){

  return(str_c(" --seqBias", # Perform sequence-specific bias correction
               " --gcBias", # Salmon will enable it to learn and correct for fragment-level GC biases in the input data
               " --posBias", # Enables modelling of position-specific fragment start distribution. This feature is considered "experimetnal", but according to Seb this has been the case since Salmon released in 2017. Recommended use anyway.
               " --validateMappings", # Performs alignment-based verification to ensure quasi-mappings give rise to a reasonable alignment before they are used for further quantification.
               " --rangeFactorizationBins 4", # The range-factorization feature allows using a data-driven likelihood factorization, which can improve quantification accuracy on certain classes of “difficult” transcripts. Recommended value is 4, as used in the range-factorisation paper: https://academic.oup.com/bioinformatics/article/33/14/i142/3953977#118769759
               " --useVBOpt", # Use the Variational Bayesian EM to optimise abundance estimates. Default option in Salmon.
               " --numBootstraps 30" # Number of bootstrap samples to generate for assessment of technical variance in the main abundance estimate produced. Useful for downstream tools (e.g. differential expression) that can use these uncertainty estimates. 30 based on Seb recommendation.
  ))

}
