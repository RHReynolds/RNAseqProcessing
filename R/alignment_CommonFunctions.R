#' Function to make a data frame of trimmed paths for STAR.
#'
#' @param fastq_dir_paths Path names in which trimmed fastq files for alignment are located.
#' @param prefix_to_sample_name Specify if samples have been given a prefix to
#'   their original name.
#' @param to_exclude_to_get_sample_name The text (or regex) that needs to be excluded
#'   from the tail end of the filename to get the sample name. If no argument
#'   provided, defaults to empty string.
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
#' @param to_exclude_to_get_sample_name The text (or regex) that needs to be excluded
#'   from the tail end of the filename to get the sample name. If no argument
#'   provided, defaults to empty string.
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
#' @param sj_dir_path Path to directory containing splice junctions files (SJ.out.tab).
#'
#' @return
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
