#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

main <- function() {

    argv <- parse_input_params()
    df <- read_tsv(argv$in_data_fp, col_types=cols())

    split_fields_list <- strsplit(df[[argv$target_column]], argv$delim)
    
    if (argv$trim_pattern != "") {
        split_fields_list <- lapply(
            split_fields_list,
            function(fields, trim_pattern) {
                gsub(trim_pattern, "", fields)
            },
            trim_pattern=argv$trim_pattern
        )
    }
    
    df[[argv$target_column]] <- lapply(
        split_fields_list,
        function(fields) {
            paste(sort(fields), collapse=argv$delim)
        }
    ) %>% unlist()
    
    write_tsv(df, path=argv$out_data_fp)
}

parse_input_params <- function() {

    parser <- arg_parser("Sort fields within selected column on IDs")

    parser <- add_argument(parser, "--in_data_fp", help="Raw data matrix file", type="character", nargs=1)
    parser <- add_argument(parser, "--out_data_fp", help="Output data matrix file", type="character", nargs=1)
    parser <- add_argument(parser, "--target_column", help="Header name of target column", type="character")
    parser <- add_argument(parser, "--delim", help="Delimitor used to distinguish fields within the target column", type="character")

    parser <- add_argument(parser, "--trim_pattern", help="Optional pattern to trim parts inside each ID (between splitters)", type="character", default="")
    
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
