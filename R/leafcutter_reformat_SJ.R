#' Convert SJ.out.tab files to leafcutter .junc file
#'
#' This function will convert SJ.out.tab files outputted by STAR alignment to
#' the .junc file format necessary for leafcutter splicing analyses.
#'
#' @param sj_dir_path Path to directory containing splice junctions files
#'   (SJ.out.tab).
#' @param output_path Path to output directory.
#'
#' @return Leafcutter formatted .junc files in the output_path.
#' @export
#'

convert_STAR_SJ_to_junc <- function(sj_dir_path, output_path){

  library(tidyverse)

  paths <- list.files(path = sj_dir_path, pattern = "_SJ.out.tab", full.names = TRUE)

  for(i in 1:length(paths)){

    sample_name <- paths[i] %>%
      str_replace("/.*/", "") %>%
      str_replace("_SJ.out.tab", "")

    cat("Loading splice junctions from:", paths[i],"\n")
    sj_out <- read_delim(file = paths[i],
                        delim = "\t",
                        col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "in_annotation", "unique_reads_junction", "multi_map_reads_junction", "max_splice_alignment_overhang"),
                        col_types = cols(chr = "c", .default = "d"))

    sj_out_leafcutter_format <-
      sj_out %>%
      dplyr::mutate(strand = ifelse(strand == 0, ".",
                                    ifelse(strand == 1, "+", "-")),
                    unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
      dplyr::select(chr, intron_start, intron_end, unknown, unique_reads_junction, strand)

    # change the stop and start into characters to avoid saving with the scientific notation
    write_delim(sj_out_leafcutter_format %>%
                  dplyr::mutate(intron_start = as.character(intron_start),
                                intron_end = as.character(end),
                                unique_reads_junction = as.character(unique_reads_junction)),
                str_c(output_path, "/", sample_name, "_SJ_leafcutter.junc"),
                delim = "\t",
                col_names = F)
  }

}




