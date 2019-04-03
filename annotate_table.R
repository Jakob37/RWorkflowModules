#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

# Procedure:
# Retrieve protein data information, for instance from UniProt
# Retrieve it in .tab format
# Run the script

main <- function() {
    
    argv <- parse_input_params()
    # source("~/src/RWorkflowModules/debug_tools.R")
    # debug_tools$use_print_argv(argv)
    
    rdf <- readr::read_tsv(argv$rdf, col_types=readr::cols(), comment = "#", na = argv$na_val)
    db_df <- readr::read_tsv(argv$db, col_types=readr::cols())

    rdf_first_protein_ids <- unlist(lapply(strsplit(rdf[[argv$rdf_name_col]], argv$rdf_name_col_splitter), function(entry) { entry[1] }))
    rdf_annotated_protein_ids <- paste0(argv$add_prefix, rdf_first_protein_ids, argv$add_suffix)
    
    match_status_message(rdf_annotated_protein_ids, db_df[[argv$db_name_col]])
    
    db_df[[argv$db_name_col]] <- as.character(db_df[[argv$db_name_col]])
    
    queries <- data.frame(Protein=rdf_annotated_protein_ids)
    queries$Protein <- as.character(queries$Protein)
    
    protein_annot_df <- queries %>% left_join(db_df, setNames(argv$db_name_col, "Protein"), keep=TRUE)
    out_df <- cbind(rdf, protein_annot_df)
    write_tsv(out_df, path=argv$out)
}

match_status_message <- function(query_proteins, db_names) {
    match_count <- length(which(query_proteins %in% db_names))
    message(
        "Match rate against database (column: \"", argv$db_name_col, "\"): ", 
        round(100 * match_count / length(protein_search_queries), 3), "% ", 
        "(", match_count, "/", length(protein_search_queries), ")")
}

parse_input_params <- function() {
    
    parser <- arg_parser("Correlation calculations")
    parser <- add_argument(parser, "--rdf", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--db", help="Tab-separated database file", type="character")
    parser <- add_argument(parser, "--out", help="Output matrix path", type="character")
    
    parser <- add_argument(parser, "--rdf_name_col", help="Name for column containing ID in rdf", type="character")
    
    parser <- add_argument(parser, "--rdf_name_col_splitter", help="Divider when multiple annotations", type="character", default=",")
    
    parser <- add_argument(parser, "--db_name_col", help="Name for column containing ID in db", type="character")
    parser <- add_argument(parser, "--db_annot_col", help="Annotation columns in db", type="character")

    parser <- add_argument(parser, "--add_prefix", help="Optional extension of input data for protein matching", type="character", default=NULL)
    parser <- add_argument(parser, "--add_suffix", help="Optional extension of input data for protein matching", type="character", default=NULL)
    parser <- add_argument(parser, "--na_val", help="NA value field", default="NA", type="character")

    argv <- parse_args(parser)
    argv
}

if (!interactive()) {
    main()
}