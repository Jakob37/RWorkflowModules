#!/usr/bin/env Rscript

library(argparser)
library(rlang)

main <- function() {
    
    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(readr))
    
    ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols())
    annot_df <- readr::read_tsv(argv$annot_fp, col_types=cols())

    col_ddf <- argv$map_col_ddf
    col_annot <- argv$map_col_annot
    
    by_vector <- c(col_ddf)
    names(by_vector) <- col_annot
    
    annotated_ddf <- ddf %>% dplyr::left_join(., annot_df, by=by_vector)
    message("Writing to: ", argv$out_fp)
    write_tsv(annotated_ddf, path = argv$out_fp)
}

parse_input_params <- function() {

    parser <- arg_parser("Append columns to design matrix using one field as a mapper")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character", nargs=1)
    parser <- add_argument(parser, "--annot_fp", help="Matrix containing additional information", type="character", nargs=1)
    parser <- add_argument(parser, "--out_fp", help="Output filepath", type="character", nargs=1)

    parser <- add_argument(parser, "--map_col_ddf", help="Name of mapping col in ddf", type="character", nargs=1)
    parser <- add_argument(parser, "--map_col_annot", help="Name of mapping col in annot", type="character", nargs=1)

    parser <- add_argument(parser, "--show_param", help="Display help output", type="bool", default=FALSE)
    parser <- add_argument(parser, "--debug_tools_path", help="Display help output", type="character", default="RWorkflowModules/debug_tools.R")
    
    argv <- parse_args(parser)
    
    if (argv$show_param) {
        source(argv$debug_tools_path)
        debug_tools$use_print_argv(argv)
    }
    
    argv
}

if (!interactive()) {
    main()
}

