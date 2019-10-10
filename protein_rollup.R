#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)
library(ProteinRollup)

# "Rscript RWorkflowModules/protein_rollup.R --ddf_fp design.tsv --sample_col sample --protein_col External.IDs --rdf_fp data.tsv --out_fp out.tsv"

main <- function() {
    
    argv <- parse_input_params()
    raw_rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols())
    
    if (argv$min_presence < 0 || argv$min_presence > 1) {
        stop("--min_presence expected to be in range 0 to 1, found: ", argv$min_presence)
    }
    
    if (!argv$one_column_mode) {
        ddf <- read_tsv(argv$ddf_fp, col_types=cols())
        rdf <- raw_rdf
        sdf <- rdf %>% dplyr::select(dplyr::one_of(ddf[[argv$sample_col]])) %>% as.matrix()
        protein_data <- rdf %>% dplyr::select(argv$protein_col) %>% unlist()
        protein_data[is.na(protein_data)] <- "NA"
    }
    else {
        sdf <- raw_rdf[, -1] %>% as.matrix()
        protein_data <- raw_rdf[, 1] %>% unlist()
        protein_data[is.na(protein_data)] <- "NA"
    }
    
    message("Performing rollup for ", nrow(raw_rdf), " peptides for ", length(unique(protein_data)), " unique protein IDs")
    
    prot_adf <- ProteinRollup::protein_rollup(
        protein_data,
        as.matrix(sdf),
        one_hit_wonders=argv$one_hit_wonders,
        min_presence=argv$min_presence,
        min_overlap=argv$min_overlap,
        debug_protein=argv$debug_protein
    )
    
    if (argv$out_protein_name != "") {
        prot_adf <- cbind(prot_adf[["Protein"]], prot_adf)
        colnames(prot_adf) <- c(argv$out_protein_name, colnames(prot_adf)[-1])
    }
    
    message(sprintf("Writing %s proteins to %s", nrow(prot_adf), argv$out_fp))
    write_tsv(prot_adf, path=argv$out_fp)
}

parse_input_params <- function() {
    
    parser <- argparser::arg_parser("Protein rollup, standalone wrapper")
    parser <- argparser::add_argument(parser, "--ddf_fp", help="Design matrix path", type="character", default=NA)
    parser <- argparser::add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character", nargs=1)
    parser <- argparser::add_argument(parser, "--out_fp", help="Output matrix path", type="character", nargs=1)
    parser <- argparser::add_argument(parser, "--sample_col", help="Design matrix sample column", type="character", default=NA)
    parser <- argparser::add_argument(parser, "--protein_col", help="Protein column in main data frame", type="character", default=NA)
    
    parser <- argparser::add_argument(parser, "--one_column_mode", help="If first column is the only annotation column, the --ddf_fp, --sample_col and --protein_col argument can be omitted", type="boolean", default=FALSE)
    
    parser <- argparser::add_argument(parser, "--min_overlap", help="Min. shared overlap between reference peptides and other", type="numeric", default=3)

    parser <- argparser::add_argument(parser, "--one_hit_wonders", help="Should proteins with a single peptide as support be included", type="boolean", default=FALSE)
    parser <- argparser::add_argument(parser, "--min_presence", help="Skip proteins where no peptide is present in at least this percentage", type="numeric", default=0)
    parser <- argparser::add_argument(parser, "--out_protein_name", help="Name of protein column in output", type="character", default="")

    parser <- argparser::add_argument(parser, "--debug_protein", help="Specify protein ID to debug", default=NA, type="character")
    parser <- argparser::add_argument(parser, "--show_warnings", help="Immediately print warnings", type="bool", default=FALSE)
    
    argv <- argparser::parse_args(parser)
    
    if (!argv$one_column_mode && (is.na(argv$sample_col) || is.na(argv$protein_col) || is.na(argv$ddf_fp))) {
        stop("If not running one_column_mode both --ddf_fp, --sample_col and --protein_col needs to be supplied")
    }
    
    if (argv$show_warnings) {
        options(warn=1)
        message("warn=1 assigned, warnings shown as they occur")
    }
    
    argv
}

if (!interactive()) {
    main()
}

