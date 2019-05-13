#!/usr/bin/env Rscript

library(argparser)

main <- function() {
    
    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(gage))
    suppressPackageStartupMessages(library(pathview))
    suppressPackageStartupMessages(library(tidyverse))
    
    if (argv$set_type == "GO") {
        # Format: Arabidopsis
        feature_set <- gage::go.gsets(species=argv$species)$go.sets
    }
    else if (argv$set_type == "KEGG") {
        # Format: ath
        feature_set <- gage::kegg.gsets(species=argv$species)$kg.sets
    }
    else {
        stop("Unknown set type: ", argv$set_type, " acceptible are 'GO' and 'KEGG'")
    }
    
    ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols())
    raw_rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols(), comment = "#", na = argv$na_val)
    rdf <- raw_rdf
    rdf$annot <- rdf %>% dplyr::select(argv$annot_col) %>% unlist()
    if (!is.na(argv$trim_reg)) {
        message("Trimming using pattern: ", argv$trim_reg)
        rdf$annot <- rdf$annot %>% gsub(argv$trim_reg, "", .)
    }
    
    # WARNING: HOW TO HANDLE REDUNDANT IDS?
    warning("Wildly removing redundant IDs with no regards for accuracy or optimal selection")
    # Could I use the ones with the most values? Then highest intensity?
    
    rdf <- rdf %>% filter(!is.na(annot))
    cond_col <- ddf %>% dplyr::select(argv$cond_col) %>% unlist() %>% as.character()
    sdf_raw <- rdf %>% dplyr::select(ddf[[argv$sample_col]]) %>% as.data.frame()
    sdf_raw$annot <- rdf$annot
    
    sdf <- reduce_dataframe_for_go(sdf_raw, ddf[[argv$sample_col]], "annot")
    # rownames(sdf) <- rdf$annot
    
    
    # rdf <- rdf %>% filter(!is.na(annot)) %>% distinct(annot, .keep_all=TRUE)
    # cond_col <- ddf %>% dplyr::select(argv$cond_col) %>% unlist() %>% as.character()
    # sdf <- rdf %>% dplyr::select(ddf[[argv$sample_col]]) %>% as.data.frame()
    # rownames(sdf) <- rdf$annot
    
    deparsed_contrasts <- deparse_contrasts(argv$contrasts, argv$cond_col, splitter="-")

    used_names <- c()
    gage_dfs <- list()    
    for (contrast_i in seq_len(nrow(deparsed_contrasts))) {
        message("Processing contrast: ", argv$contrasts[[contrast_i]])
        contrast_levels <- deparsed_contrasts[contrast_i, ]
        
        ref_inds <- which(cond_col == contrast_levels[1])
        samp_inds <- which(cond_col == contrast_levels[2])
        
        if (length(ref_inds) == 0 || length(samp_inds) == 0) {
            message("At least one level didn't have samples, skipping. ref_inds entries: ", length(ref_inds), " samp_inds entries: ", length(samp_inds))
            next
        }
        
        gage_out <- gage::gage(
            sdf,
            gsets=feature_set,
            ref=ref_inds,
            samp=samp_inds,
            compare="unpaired",
            same.dir=TRUE
        )

        message("Processing contrast_i: ", contrast_i)
        used_names <- c(used_names, argv$contrast_names[contrast_i])        
        gage_df <- gage_out %>% as.data.frame() %>% rownames_to_column() %>% rowid_to_column()
        gage_dfs[[argv$contrast_names[contrast_i]]] <- gage_df
    }
    names(gage_dfs) <- used_names
    combined_dfs <- do.call("cbind", gage_dfs)
    
    write_tsv(combined_dfs, path=argv$enrichment_out_fp)
}

reduce_dataframe_for_go <- function(df, samples, reduce_col) {
    
    unique_levels <- unique(df[[reduce_col]])
    unique_lines <- list()
    
    for (unique_level in unique_levels) {
        unique_lines[[unique_level]] <- df %>% 
            filter(annot == unique_level) %>% 
            dplyr::select(samples) %>% 
            mutate(select_col=rowSums(., na.rm=TRUE)) %>% 
            filter(select_col==max(select_col)) %>%
            head(1) %>%
            dplyr::select(-select_col)
    }
    do.call("rbind", unique_lines)
}

deparse_contrasts <- function(contrasts, contrast_base, splitter="-") {
    splits <- do.call("rbind", strsplit(contrasts, splitter))
    trimmed_splits <- gsub(contrast_base, "", splits)
    trimmed_splits
}

parse_input_params <- function() {
    
    parser <- arg_parser("Enrichment module")
    parser <- add_argument(parser, "--rdf_fp", help="Data matrix path", type="character")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
    parser <- add_argument(parser, "--sample_col", help="Sample column header in design matrix", type="character")
    parser <- add_argument(parser, "--cond_col", help="Condition column header in design matrix", type="character")
    parser <- add_argument(parser, "--annot_col", help="Annot. col in data matrix", type="character")
    parser <- add_argument(parser, "--enrichment_out_fp", help="Output path for annotated matrix", type="character")
    
    parser <- add_argument(parser, "--contrasts", help="Condition contrasts for which to perform enrichment, two entries, first is reference level, second studied condition", type="character", nargs=Inf)
    parser <- add_argument(parser, "--contrast_names", help="Names of provided contrasts", type="character", nargs=Inf)
    
    parser <- add_argument(parser, "--species", help="Target organism in GAGE database", type="character")
    
    parser <- add_argument(parser, "--set_type", help="Set type, either GO (default) or KEGG", type="character", default="GO")
    parser <- add_argument(parser, "--trim_reg", help="Trim parts of ID matching pattern", type="character", default=NA)
    parser <- add_argument(parser, "--na_val", help="NA value in input", type="character", default="NA")
    
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

