#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

message("Executing")

# Procedure:
# Retrieve protein data information, for instance from UniProt
# Retrieve it in .tab format
# Run the script

main <- function() {
    
    argv <- parse_input_params()

    rdf <- readr::read_tsv(argv$rdf, col_types=readr::cols(), comment = "#", na = argv$na_val)
    db_df <- readr::read_tsv(argv$db, col_types=readr::cols())

    rdf_annotated_protein_ids <- str_split(rdf[[argv$rdf_name_col]], argv$rdf_name_col_splitter, simplify=TRUE)[, 1]
    db_df[[argv$db_name_col]] <- as.character(db_df[[argv$db_name_col]])
    match_status_message(rdf_annotated_protein_ids, db_df[[argv$db_name_col]], argv$db_name_col)
    
    queries <- data.frame(Protein=rdf_annotated_protein_ids)
    queries$Protein <- as.character(queries$Protein)
    
    protein_annot_df <- queries %>% left_join(db_df, setNames(argv$db_name_col, "Protein"), keep=TRUE)
    out_df <- cbind(rdf, protein_annot_df)
    write_tsv(out_df, path=argv$out)
}

match_status_message <- function(query_proteins, db_names, db_col_name) {
    match_count <- length(which(query_proteins %in% db_names))
    message(
        "Match rate against database (column: \"", db_col_name, "\"): ", 
        round(100 * match_count / length(query_proteins), 3), "% ", 
        "(", match_count, "/", length(query_proteins), ")")
}

parse_input_params <- function() {
    
    parser <- arg_parser("Correlation calculations")
    parser <- add_argument(parser, "--rdf", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--db", help="Tab-separated database file", type="character")
    parser <- add_argument(parser, "--out", help="Output matrix path", type="character")
    
    parser <- add_argument(parser, "--rdf_name_col", help="Name for column containing ID in rdf", type="character")
    parser <- add_argument(parser, "--rdf_name_col_splitter", help="Divider when multiple annotations", type="character", default=",")
    parser <- add_argument(parser, "--db_name_col", help="Name for column containing ID in db", type="character")
    # parser <- add_argument(parser, "--db_annot_col", help="Annotation columns in db", type="character")

    parser <- add_argument(parser, "--na_val", help="NA value field", default="NA", type="character")

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