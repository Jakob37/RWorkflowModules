#!/usr/bin/env Rscript

library(argparser)

message("Requires NormalyzerDE to run (github.com/ComputationalProteomics/NormalyzerDE)")

main <- function() {
    
    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(readr))
    
    ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols())
    
    if (!argv$sample_col %in% colnames(ddf)) {
        stop("Sample column must be present in design matrix header, looking for column: ", argv$sample_col,
             "\nFound: ", paste(colnames(ddf), collapse=", "))
    }
    
    rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols(), comment = "#", na = argv$na_val)
    sdf <- rdf %>% dplyr::select(ddf[[argv$sample_col]]) %>% as.matrix()
    adf <- rdf %>% dplyr::select(-one_of(ddf[[argv$sample_col]]))
    
    message("Loaded raw data with dimensions: ", paste(dim(rdf), collapse=", "))

    if (argv$normalization == "raw") {
        norm_method <- NULL
    }
    else if (argv$normalization == "loess") {
        norm_method <- NormalyzerDE::performCyclicLoessNormalization
    }
    else if (argv$normalization == "median") {
        norm_method <- NormalyzerDE::medianNormalization
    }
    else if (argv$normalization == "quantile") {
        norm_method <- NormalyzerDE::performQuantileNormalization
    }
    else if (argv$normalization == "log2") {
        norm_method <- log2
    }
    else if (argv$normalization == "vsn") { 
        if (argv$do_log2) {
            stop("VSN normalization does not expect log2-transformed input, incompatible options!")
        }
        norm_method <- NormalyzerDE::performVSNNormalization
    }
    else {
        stop("Unknown normalization: ", argv$normalization, ", acceptible are: none, loess, median, vsn")
    }
    
    generate_normalized_data(
        ddf, 
        sdf, 
        adf, 
        argv$sample_col, 
        argv$out_fp, 
        norm_method=norm_method, 
        do_log2=argv$do_log2, 
        do_zscale=argv$do_zscale,
        do_groupwise=argv$do_group_norm,
        group_col=argv$group_col
    )
}

scale_genes_to_z_dist <- function(sdf) {
    
    # Equivalent in this case    
    # z_scale <- function(row) {
    #     (row - mean(row, na.rm=TRUE)) / sd(row, na.rm=TRUE)
    # } 
    # scaled_mat <- apply(sdf, 1, z_scale)
    
    scaled_mat <- t(scale(t(sdf)))
    scaled_mat
}

generate_normalized_data <- function(ddf, sdf, adf, sample_col, out_fp, norm_method=NULL, do_log2=FALSE, do_zscale=FALSE, do_groupwise=FALSE, group_col=NULL) {
    
    if (do_log2) {
        sdf <- log2(sdf)
    }

    if (!is.null(norm_method)) {
        if (!do_groupwise) {
            norm_sdf <- norm_method(sdf)
        }
        else {
            if (is.null(group_col)) {
                stop("Using 'do_groupwise' option requires that you specify the 'group_col' argument")
            }
            norm_sdf <- groupwise_normalized_data(ddf, sdf, group_col, norm_method)
        }
    }
    else {
        norm_sdf <- sdf
    }
    
    if (do_zscale) {
        norm_sdf <- scale_genes_to_z_dist(norm_sdf)
    }
    
    norm_rdf <- cbind(adf, norm_sdf)
    message("Writing ", nrow(norm_rdf), " rows to ", out_fp)
    readr::write_tsv(norm_rdf, path=out_fp)
}

groupwise_normalized_data <- function(ddf, sdf, group_col, norm_method) {
    
    groups <- ddf[[group_col]]
    unique_groups <- sort(unique(groups))

    group_indices <- lapply(unique_groups, function(unique_val, groups) {
        which(unique_val == groups)
    }, groups=groups)
    names(group_indices) <- unique_groups
    
    sample_subsets <- lapply(group_indices, function(indices, sdf) { sdf[, indices] }, sdf=sdf)
    normalized_subsets <- lapply(sample_subsets, function(sample_subset, norm) { norm(sample_subset) }, norm=norm_method)
    
    group_normalized_mat <- do.call("cbind", normalized_subsets)
    orig_order_normalized_mat <- group_normalized_mat[, unlist(group_indices)]
    
    orig_order_normalized_mat
}

parse_input_params <- function() {
    parser <- arg_parser("Normalization wrapper for batch analysis")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
    parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
    parser <- add_argument(parser, "--na_val", help="NA value fields", type="character", default="NA")
    
    parser <- add_argument(parser, "--do_group_norm", help="Perform group-wise normalization (requires ddf group column)", default=FALSE)
    parser <- add_argument(parser, "--group_col", help="Group column")
        
    parser <- add_argument(parser, "--do_zscale", help="Do z-scaling before normalization", type="logical", default=FALSE)
    parser <- add_argument(parser, "--do_log2", help="Do log-transform", type="logical", default=TRUE)
    
    parser <- add_argument(parser, "--out_fp", help="Output matrix path", type="character", default="")
    parser <- add_argument(parser, "--normalization", help="Normalization type (raw, log2, loess, median, vsn)", type="character")

    parser <- add_argument(parser, "--show_pars", help="Display help output", type="bool", default=FALSE)
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

