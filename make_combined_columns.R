#!/usr/bin/env Rscript

library(argparser)
library(rlang)
suppressPackageStartupMessages(library(tidyverse))

main <- function() {
    
    argv <- parse_input_params()
    ddf <- read_tsv(argv$in_fp, col_types=cols())
    message(sprintf("Loaded raw matrix with dimensions: %s from %s", paste(dim(ddf), collapse=","), argv$in_fp))
    
    add_cols <- do.call(
        "cbind", 
        lapply(
            argv$base_cols, 
            function(base_col, add_col, ddf) {
                paste(ddf[[base_col]], ddf[[add_col]], sep="_")
            }, 
            add_col=argv$add_col,
            ddf=ddf
        )
    )
    
    if (!is.na(argv$add_name))          add_col_name <- argv$add_name
    else                                add_col_name <- argv$add_col

    if (!is.na(argv$base_col_names))    base_col_names <- argv$base_col_names
    else                                base_col_names <- argv$base_cols
        
    colnames(add_cols) <- paste(base_col_names, add_col_name, sep="_")
    ddf_added <- cbind(ddf, add_cols)
    message(sprintf("Writing updated matrix with dimensions: %s to %s", paste(dim(ddf_added), collapse=","), argv$out_fp))
    write_tsv(ddf_added, path=argv$out_fp)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Sort fields within selected column on IDs")
    
    parser <- add_argument(parser, "--in_fp", help="Raw data matrix file", type="character", nargs=1)
    parser <- add_argument(parser, "--out_fp", help="Output data matrix file", type="character", nargs=1)

    parser <- add_argument(parser, "--base_cols", help="Column with sample name matching expression data", type="character", nargs=Inf)
    parser <- add_argument(parser, "--add_col", help="Column with the base sample number", type="character", nargs=1)

    parser <- add_argument(parser, "--base_col_names", help="Base col names (instead of --base_cols in column names), must be same length as --base_cols", nargs=Inf, default=NA)
    parser <- add_argument(parser, "--add_name", help="Replacement name used in headers (instead of the column name from --add_col)", nargs=1, default=NA)
    
    parser <- add_argument(parser, "--show_param", help="Show input parameters, for debug", type="bool", default=FALSE)
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
