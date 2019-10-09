#!/usr/bin/env Rscript

library(argparser)
library(rlang)

main <- function() {
    
    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(readr))

    raw_rdf <- read_tsv(argv$rdf_fp, col_types=cols())
    
    add_dfs <- lapply(strsplit(argv$annot_fps, " "), function(annot_fp) {
        read_tsv(annot_fp, col_types=cols())
    })

    add_df <- do.call("cbind", add_dfs)

    if (argv$ddf_fp == "") {
        out_df <- cbind(add_df, raw_rdf)
    } else {
        ddf <- read_tsv(argv$ddf, col_types=cols())
        sdf <- raw_rdf %>% dplyr::select(ddf[[argv$sample_col]])
        adf <- raw_rdf %>% dplyr::select(-one_of(ddf[[argv$sample_col]]))
        out_df <- cbind(adf, add_df, sdf)
    }
    
    message("Writing to: ", argv$out_fp)
    write_tsv(out_df, path = argv$out_fp)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Insert annotation columns into rdf data frame")
    parser <- add_argument(parser, "--rdf_fp", help="Raw data matrix", type="character")
    parser <- add_argument(parser, "--annot_fps", help="Space-separated list of annotation files to insert", type="character", nargs=Inf)
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix, can be used to insert columns after initial annotations", type="character", default="")
    parser <- add_argument(parser, "--sample_col", help="Specify what column in design matrix contains samples", type="character")
    parser <- add_argument(parser, "--out_fp", help="Output filepath", type="character")
    
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

