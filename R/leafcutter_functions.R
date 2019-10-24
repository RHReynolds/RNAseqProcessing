#' Convert SJ.out.tab files to leafcutter .junc file and remove ENCODE blacklist
#' regions
#'
#' This function will convert SJ.out.tab files outputted by STAR alignment to
#' the .junc file format necessary for leafcutter splicing analyses. THis
#' function will also remove any junctions that overlap with ENCODE blacklist
#' regions (https://www.nature.com/articles/s41598-019-45839-z)
#'
#' @param sj_dir_path Path to directory containing splice junctions files
#'   (SJ.out.tab).
#' @param output_path Path to output directory where .junc files will be stored.
#' @param filter_out_blacklist_regions To filter out ENCODE blacklist regions,
#'   set to TRUE. Default = FALSE.
#' @param path_to_ENCODE_blacklist If filter_out_blacklist_regions = TRUE, must
#'   provide path to a bed file with ENCODE blacklist regions
#'   (https://github.com/Boyle-Lab/Blacklist/tree/master/lists).
#'
#' @return Leafcutter-formatted .junc files in the output_path, in addition to a
#'   list of generated .junc files in a separate .txt file.
#' @export
#'

convert_STAR_SJ_to_junc <- function(sj_dir_path, output_path, filter_out_blacklist_regions=FALSE, path_to_ENCODE_blacklist = NULL){

  library(tidyverse)

  paths <- list.files(path = sj_dir_path, pattern = "_SJ.out.tab", full.names = TRUE)

  if(filter_out_blacklist_regions == TRUE){

    library(GenomicRanges)
    library(rtracklayer)

    # Load encode blacklist (https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
    ENCODE_blacklist_hg38 <- rtracklayer::import(path_to_ENCODE_blacklist)

  }

  for(i in 1:length(paths)){

    sample_name <- paths[i] %>%
      str_replace("/.*/", "") %>%
      str_replace("_SJ.out.tab", "")

    cat("Loading splice junctions from:", paths[i],"\n")

    sj_out <- read_delim(file = paths[i],
                         delim = "\t",
                         col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "in_annotation", "unique_reads_junction", "multi_map_reads_junction", "max_splice_alignment_overhang"),
                         col_types = cols(chr = "c", .default = "d"))

    if(filter_out_blacklist_regions == TRUE){

      sj_out <- sj_out %>%
        dplyr::mutate(chr = str_c("chr", chr),
                      strand = ifelse(strand == 0, "*",
                                      ifelse(strand == 1, "+", "-")),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        makeGRangesFromDataFrame(.,
                                 keep.extra.columns = TRUE,
                                 seqnames.field = "chr",
                                 start.field = "intron_start",
                                 end.field = "intron_end",
                                 ignore.strand = FALSE)

      # Remove junctions that overlap with ENCODE blacklist regions
      overlapped_junctions <- GenomicRanges::findOverlaps(query = ENCODE_blacklist_hg38,
                                                          subject = sj_out,
                                                          ignore.strand = F)

      indexes <- subjectHits(overlapped_junctions)
      sj_out <- sj_out[-indexes, ] %>%
        as.data.frame() %>%
        dplyr::rename(chr = seqnames,
                      intron_start = start,
                      intron_end = end)

      # Convert to leafcutter format
      sj_out_leafcutter_format <-
        sj_out %>%
        dplyr::mutate(chr = str_remove(chr, "chr"),
                      strand = str_replace(strand, "\\*", "."),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        dplyr::select(chr, intron_start, intron_end, unknown, unique_reads_junction, strand)

    } else{

      sj_out_leafcutter_format <-
        sj_out %>%
        dplyr::mutate(strand = ifelse(strand == 0, ".",
                                      ifelse(strand == 1, "+", "-")),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        dplyr::select(chr, intron_start, intron_end, unknown, unique_reads_junction, strand)

    }

    # change the stop and start into characters to avoid saving with the scientific notation
    write_delim(sj_out_leafcutter_format %>%
                  dplyr::mutate(intron_start = as.integer(intron_start),
                                intron_end = as.integer(intron_end),
                                unique_reads_junction = as.integer(unique_reads_junction)),
                path = str_c(output_path, "/", sample_name, "_SJ_leafcutter.junc"),
                delim = "\t",
                col_names = F)
  }

  # write a .txt file with each
  junc_df <- tibble(junc_file_name = list.files(path = output_path, pattern = "_SJ_leafcutter.junc", full.names = TRUE))

  write_delim(junc_df,
              path = str_c(output_path, "/list_juncfiles.txt"),
              delim = "\t",
              col_names = F)

}



#' Create LeafCutter group files for multiple pairwise comparisons.
#'
#' LeafCutter's differential splicing analyses currently only support pairwise
#' comparisons. For each pairwise comparison it requires a group file to specify
#' which samples belong to which group. Thus, if a grouping variable contains > 2
#' groups, multiple pairwise comparisons must be made. This function will
#' identify the comparisons, based on an inputted grouping column, and will
#' output separate group .txt files for each group comparison combination.
#'
#' @param df Dataframe that should have a minimum of two columns in the
#'   following order: (i) sample name (filename in the leafcutter
#'   '_perind_numers.count.gz' file) and (ii) grouping variable. Additional
#'   columns will be used as confounders.
#' @param group_column_name Name of column with grouping variable. Must be
#'   entered in quotation marks "".
#' @param output_path Output path where _group_file.txt files will be stored for
#'   later LeafCutter differential splicing analyses. Files will be named using
#'   group names for the comparison in question.
#'
#' @return Named _group_file.txt files for LeafCutter differential splicing
#'   analyses. Files will be named using the groups included in the pairwise
#'   comparison.
#' @export
#'

create_group_files_multi_pairwisecomp <- function(df, group_column_name, output_path){

  library(dplyr)
  library(lazyeval)

  # Create vector of groups
  groups <- df %>%
    .[[group_column_name]] %>%
    unique()

  # Create dataframe of comparisons
  comparisons <-
    combn(x = groups,
          m = 2) %>%
    t() %>%
    as_tibble()

  # Loop through comparisons and for each write out a group file
  for(i in 1:nrow(comparisons)){

    # Filter for comparisons
    comparison <- comparisons[i,] %>%
      as_vector()

    # Create filter expression using column name and comparison vector
    filter_criteria <- interp(~y %in% x, .values = list(y = as.name(group_column_name), x = comparison))

    select_comparison <- df %>%
      dplyr::filter_(filter_criteria)

    write_delim(x = select_comparison,
                path = str_c(output_path, "/",
                             comparison[1], "_vs_", comparison[2], "_",
                             "group_file.txt"),
                delim = "\t",
                col_names = F)

  }

}

#' Run Leafviz app
#'
#' Function that allows leafviz shiny app to run when called from inside
#' RStudio.
#'
#' @param leafviz_dir Path to the directory containing the (i) the
#'   run_leafviz.R, (ii) the server.R and (iii) ui.R scripts.
#' @param results_filename Name of the file (including the path to it) that
#'   should be loaded by leafviz.
#'
#' @return The Leafviz app will open in a separate browser and load the results
#'   of a differential splicing analysis.
#' @export
#'

run_leafviz <- function(leafviz_dir, results_filename){

  # Load libraries
  library(shiny, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(DT, quietly = TRUE)
  library(leafcutter,quietly = TRUE)
  library(reshape2, quietly = TRUE)
  library(gridExtra, quietly = TRUE)
  library(intervals, quietly = TRUE) # needed for pretty strand arrow placement
  library(ggplot2, quietly = TRUE)
  library(foreach, quietly = TRUE)
  library(shinycssloaders, quietly = TRUE)

  # Set working directory
  setwd(leafviz_dir)

  # Load results file
  if(length(results_filename) == 0){
    print("No results found!")
  }

  print(paste0("Loading results from ",results_filename))
  load(results_filename, envir=.GlobalEnv)

  # Start shiny app
  shiny::runApp( launch.browser=TRUE )


}

