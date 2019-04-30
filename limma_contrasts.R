#!/usr/bin/env Rscript

library(argparser)

main <- function() {
    
    argv <- parse_input_params()
    
    if (!is.na(argv$contrast_names) && length(argv$contrasts) != length(argv$contrast_names)) {
        stop("Argument '--contrast_names' must either be NA or of same length as '--contrasts', found: ",
        paste(argv$contrasts, collapse=", "), " and ", paste(argv$contrast_names))
    }
    
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(limma))
    
    raw_ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols())
    
    if (!is.na(argv$omit_samples) && length(argv$omit_samples) > 0) {
        ddf <- raw_ddf %>% as.data.frame() %>% dplyr::filter(!(UQ(as.name(argv$sample_col)) %in% argv$omit_samples))
    } else {
        ddf <- raw_ddf
    }
    
    rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols(), comment = "#", na = argv$na_val)
    raw_sdf <- rdf %>% dplyr::select(raw_ddf[[argv$sample_col]]) %>% as.matrix()
    sdf <- rdf %>% dplyr::select(ddf[[argv$sample_col]]) %>% as.matrix()
    adf <- rdf %>% dplyr::select(-one_of(raw_ddf[[argv$sample_col]]))
    
    message("Loaded raw data with dimensions: ", paste(dim(rdf), collapse=", "))
    
    named_contrast_vect <- setup_contrasts(ddf, argv$cond_col_name, argv$contrasts, argv$omit_absent_contrasts, argv$contrast_names)
    
    # Prepare R equation objects
    combined_limma_tables <- calculate_combined_limma_tables(argv$model, ddf, sdf, unname(named_contrast_vect), contrast_names=names(named_contrast_vect))
    
    # Output
    full_stat_rdf <- cbind(adf, combined_limma_tables, raw_sdf)
    
    message("Writing data with dimensions: ", paste(dim(full_stat_rdf), collapse=", "))
    readr::write_tsv(full_stat_rdf, path=argv$out_fp)
}

setup_contrasts <- function(ddf, cond_col_name, contrasts, omit_absent_contrasts, contrast_names=NULL) {
    
    valid_contrasts <- check_valid_contrasts(contrasts, cond_col_name, ddf)
    contrast_vect <- contrasts[valid_contrasts]

    if (length(which(!valid_contrasts)) > 0) {
        if (!omit_absent_contrasts) {
            stop("Contrasts missing, stopping. If you want to proceed - assign '--omit_absent_contrasts'")
        }
        else {
            warning(
                "Contrasts omitted: ", paste(contrasts[!valid_contrasts], collapse=", "), 
                " with names: ", paste(contrast_names[!valid_contrasts], collapse=", ")
            )
        }
    }

    contrast_vect <- contrasts[valid_contrasts]
    if (!is.null(contrast_names)) {
        names(contrast_vect) <- contrast_names[valid_contrasts]
    }
    
    contrast_vect
}

check_valid_contrasts <- function(contrasts, cond_col_name, ddf, splitter="-") {
    
    valid_contrasts <- lapply(
        strsplit(contrasts, splitter), 
        function(pair, ddf) { 
            all(gsub(cond_col_name, "", pair) %in% ddf[[cond_col_name]]) 
        }, ddf=ddf) %>% unlist()
    
    valid_contrasts
}

calculate_combined_limma_tables <- function(model_string, ddf, sdf, contrasts, contrast_names=NA) {
    
    model <- as.formula(model_string)
    model_design <- model.matrix(model, data=ddf)
    
    # Calculate Limma table (pairwise)
    fit <- limma::lmFit(sdf, model_design)
    contrast.matrix <- limma::makeContrasts(contrasts=contrasts, levels=model_design)
    
    fit_contrasts <- contrasts.fit(fit, contrast.matrix)
    fit_bayes <- limma::eBayes(fit_contrasts)
    
    limma_tables <- lapply(
        seq_len(length(colnames(fit_bayes$coefficients))), 
        function(coef) { 
            my_tbl <- topTable(fit_bayes, coef=coef, number=Inf, sort.by="none")
            my_tbl$absLogFC <- abs(my_tbl$logFC)
            my_tbl
        })

    if (all(is.na(contrast_names))) {
        names(limma_tables) <- colnames(fit_bayes$coefficients)
    }
    else {
        names(limma_tables) <- contrast_names
    }
    combined_limma_tables <- do.call("cbind", limma_tables)
    combined_limma_tables <- cbind(row_nbr=seq_len(nrow(combined_limma_tables)), combined_limma_tables)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Statistical comparison wrapper")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
    parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--out_fp", help="Output matrix path", type="character")
    parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
    parser <- add_argument(parser, "--na_val", help="NA value fields", type="character", default="NA")
    parser <- add_argument(parser, "--omit_samples", help="List of samples to omit", type="character", nargs=Inf)
    
    parser <- add_argument(parser, "--model", help="Statistical model to use", type="character")
    parser <- add_argument(parser, "--contrasts", help="Contrasts to perform", type="character", nargs = Inf)
    parser <- add_argument(parser, "--contrast_names", help="Names for contrasts, same length as --contrasts argument", type="character", nargs = Inf, default=NA)
    
    parser <- add_argument(parser, "--cond_col_name", help="Column in design matrix containing contrast levels, used to verify contrasts if specified", type="character")
    parser <- add_argument(parser, "--omit_absent_contrasts", help="If set and contrasts are missing - continue processing", type="logical", default=FALSE)
    
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
