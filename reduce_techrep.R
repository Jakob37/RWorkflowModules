#!/usr/bin/env Rscript

library(argparser)

main <- function() {

    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(tidyverse))
    
    ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols()) %>% as.data.frame()
    rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols(), comment="#", na=argv$na_val)
    
    if (!argv$sample_col %in% colnames(ddf)) {
        stop("Sample column must be present in design matrix header, looking for column: ", argv$sample_col,
             "\nFound: ", paste(colnames(ddf), collapse=", "))
    }
    
    sdf <- rdf %>% dplyr::select(as.character(ddf[[argv$sample_col]])) %>% as.matrix()
    
    if (typeof(sdf) == "character") {
        stop("Data part loaded as character format, likely indicates incorrect ",
             " sample selection from data frame, incorrect NA value parsing or invalid field content")
    }
    
    adf <- rdf %>% dplyr::select(-one_of(as.character(ddf[[argv$sample_col]]))) %>% as.data.frame()
    
    red_mats <- reduce_technical_replicates_for_matrices(designMat = ddf, dataMat = sdf, techRepGroups=ddf[[argv$techrep_col]])
    
    if (!is.na(argv$techred_rdf_fp)) {
        message("Writing ", nrow(red_mats$data), " to ", argv$techred_rdf_fp)
        write_tsv(cbind(adf, data.frame(red_mats$data)), path=argv$techred_rdf_fp)
    }
    
    if (!is.na(argv$techred_ddf_fp)) {
        message("Writing ", nrow(red_mats$design), " to ", argv$techred_ddf_fp)
        write_tsv(red_mats$design, path=argv$techred_ddf_fp)
    }
}

reduce_technical_replicates_for_matrices <- function(dataMat, designMat, techRepGroups) {
    
    uniqueGroups <- unique(techRepGroups)
    indices <- lapply(uniqueGroups, function(i) { which(techRepGroups %in% i) })
    
    collData <- lapply(
        indices, 
        function(inds) { rowMeans(dataMat[, inds, drop=FALSE], na.rm = TRUE) })
    
    collDataMat <- as.matrix(data.frame(collData))
    colnames(collDataMat) <- uniqueGroups
    
    first_indices <- sapply(indices, function(ind_group){ind_group[1]})
    collDesignMat <- designMat[first_indices, ]
    collDesignMat$sample <- uniqueGroups
    
    list(data=collDataMat, design=collDesignMat)
}

parse_input_params <- function() {
    
    message("Requires NormalyzerDE to run (github.com/ComputationalProteomics/NormalyzerDE)")
    
    parser <- arg_parser("Normalization wrapper for batch analysis")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
    parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
    parser <- add_argument(parser, "--techrep_col", help="Design matrix techrep column", type="character")
    parser <- add_argument(parser, "--na_val", help="NA value fields", type="character", default="NA")
    
    parser <- add_argument(parser, "--techred_ddf_fp", help="Techred design matrix path", type="character", default=NA)
    parser <- add_argument(parser, "--techred_rdf_fp", help="Techred data matrix path", type="character", default=NA)
    
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

