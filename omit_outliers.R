#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

main <- function() {
    
    argv <- parse_input_params()
    rdf <- read_tsv(argv$in_data_fp, col_types=cols())
    ddf <- read_tsv(argv$in_design_fp, col_types=cols())
    all_samples <- ddf[[argv$sample_col]]
    
    if (!all(argv$outliers %in% all_samples)) {
        stop("Not all outliers present in raw data, outliers: ", paste(argv$outliers, collapse=", "), " raw data headers: ", paste(colnames(rdf), collapse=", "))
    }
    
    sdf <- rdf %>% dplyr::select(all_samples)
    adf <- rdf %>% dplyr::select(-one_of(all_samples))
    
    sdf_noout <- sdf %>% dplyr::select(-one_of(argv$outliers))
    ddf_noout <- ddf %>% filter(!(UQ(as.name(argv$sample_col)) %in% argv$outliers))
    
    rdf_noout <- cbind(adf, sdf_noout)
    
    write_tsv(rdf_noout, path=argv$out_data_fp)
    write_tsv(ddf_noout, path=argv$out_design_fp)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Sort fields within selected column on IDs")
    
    parser <- add_argument(parser, "--in_data_fp", help="Raw data matrix file", type="character")
    parser <- add_argument(parser, "--out_data_fp", help="Output data matrix file", type="character")
    parser <- add_argument(parser, "--in_design_fp", help="Raw design matrix file", type="character")
    parser <- add_argument(parser, "--out_design_fp", help="Output design matrix file", type="character")
    
    parser <- add_argument(parser, "--sample_col", help="Column in design matrix with sample names", type="character")
    parser <- add_argument(parser, "--outliers", help="Outliers to remove", type="character", nargs=Inf)
    
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
