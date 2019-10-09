#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

main <- function() {
    
    argv <- parse_input_params()
    my_df <- read_tsv(argv$in_fp, col_types=cols())

    if (!(argv$target_column %in% colnames(my_df))) {
        stop("--target_column (", argv$target_column, ") not found among colnames: ", paste(colnames(my_df), collapse=", "))
    }
    
    # Vector of booleans for each set of protein IDs
    screen_grep_pattern <- paste(argv$target_patterns, collapse="|")
    pattern_matches_list <- lapply(
        strsplit(my_df[[argv$target_column]], argv$delim), 
        function(ids, screen_pat) { 
            grepl(screen_pat, ids) 
        }, 
        screen_pat=screen_grep_pattern
    )
    
    any_match <- lapply(pattern_matches_list, function(patterns) { any(patterns) }) %>% unlist()
    all_match <- lapply(pattern_matches_list, function(patterns) { all(patterns) }) %>% unlist()
    message("Total entries: ", length(any_match), " All screen match: ", table(all_match)["TRUE"], " At least one match: ", table(any_match)["TRUE"])

    if (!argv$screen_any_match) {
        filtered_df <- my_df[!all_match, ]
    }
    else {
        filtered_df <- my_df[!any_match, ]
    }
    
    message("Proceeding with --screen_any_match set to: ", argv$screen_any_match, " retaining ", nrow(filtered_df), " entries")
    write_tsv(filtered_df, path=argv$out_fp)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Sort fields within selected column on IDs")
    
    parser <- add_argument(parser, "--in_fp", help="Raw data matrix file", type="character")
    parser <- add_argument(parser, "--out_fp", help="Output data matrix file", type="character")
    parser <- add_argument(parser, "--target_column", help="Header name of target column", type="character")
    parser <- add_argument(parser, "--target_patterns", help="Comma delimited patterns to screen for", type="character", nargs=Inf)
    parser <- add_argument(parser, "--delim", help="Delimitor used to distinguish fields within the target column", type="character")
    parser <- add_argument(
        parser, "--screen_any_match", 
        help="When multiple IDs, do all IDs need to be 'unclean' for the feature to be removed, or just one", type="bool", default=FALSE)
    
    parser <- add_argument(parser, "--show_pars", help="Show input parameters, for debug", type="bool", default=FALSE)
    parser <- add_argument(parser, "--debug_tools_path", help="Display help output", type="character", default="RWorkflowModules/debug_tools.R")
    
    argv <- parse_args(parser)
    
    if (argv$show_pars) {
        source(argv$debug_tools_path)
        debug_tools$use_print_argv(argv)
    }
    
    argv
}

if (!interactive()) {
    main()
}
