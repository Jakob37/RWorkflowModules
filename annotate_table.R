#!/usr/bin/env Rscript

library(argparser)
suppressPackageStartupMessages(library(tidyverse))

message("Starting annotation script...")

# Procedure:
# Retrieve protein data information, for instance from UniProt
# Retrieve it in .tab format
# Run the script

main <- function() {
    
    argv <- parse_input_params()

    rdf <- readr::read_tsv(argv$rdf, col_types=readr::cols(), comment = "#", na = argv$na_val)
    db_df <- readr::read_tsv(argv$db, col_types=readr::cols()) %>% distinct(UQ(as.name(argv$db_name_col)), .keep_all=TRUE)
    db_df[[argv$db_name_col]] <- as.character(db_df[[argv$db_name_col]])  # Make sure database column is in string format
    
    if (argv$clean_up_descr_type == "TAIR") {
        if (is.na(argv$clean_up_descr_colname)) {
            message("No clean_up_colname assigned to clean up, skipping...")
        }
        else {
            message("Cleaning up annotations - might not be well adjusted for data! Trimming everything beyond inline-pipe character, and after 'OS='")
            db_df$Description <- db_df$Description %>% gsub("^\\| ", "", .) %>% gsub(" \\|.*", "", .) %>% gsub("OS=.*", "", .)
        }
    }
    else if (!is.na(argv$clean_up_descr_type)) {
        message("Only TAIR supported as clean_up_annot type")
    }
    
    rdf_annotated_protein_ids_list <- str_split(rdf[[argv$rdf_name_col]], argv$rdf_name_col_splitter)

    rdf_proteins_long <- tibble(id=seq_along(rdf_annotated_protein_ids_list), Protein=rdf_annotated_protein_ids_list) %>% unnest() %>% data.frame()
    protein_annot_df <- rdf_proteins_long %>% left_join(db_df, setNames(argv$db_name_col, "Protein"), keep=TRUE)

    # Retrive unique matches per line, where multiple distinct are column-wise concatenated delimited by a comma    
    out <- protein_annot_df %>% 
        group_by(id) %>% 
        group_map(
            function(df, group_info) { 
                df[, -1] %>% 
                    filter(!is.na(ProteinID)) %>% 
                    distinct(ProteinID, .keep_all = TRUE) %>% 
                    apply(2, function(col) { paste(col, collapse=",") }) %>% 
                    t() %>% 
                    data.frame() 
            })
    protein_annot_df <- do.call("rbind", out)
    protein_annot_df$ProteinIDFirst <- protein_annot_df$ProteinID %>% as.character() %>% strsplit(",") %>% lapply(function(entries) { entries[1] }) %>% unlist()
    
    message("Annotation count distribution - Number of entries with corresponding number matching proteins")
    print(protein_annot_df$ProteinID %>% unlist() %>% unname() %>% as.character() %>% strsplit(",") %>% lapply(., function(entries) { length(entries) }) %>% unlist() %>% table())
    
    message("Writing entries to file")
    out_df <- cbind(rdf, protein_annot_df)
    write_tsv(out_df, path=argv$out)
}

parse_input_params <- function() {
    
    parser <- arg_parser("Correlation calculations")
    parser <- add_argument(parser, "--rdf", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--db", help="Tab-separated database file", type="character")
    parser <- add_argument(parser, "--out", help="Output matrix path", type="character")
    
    parser <- add_argument(parser, "--rdf_name_col", help="Name for column containing ID in rdf", type="character")
    parser <- add_argument(parser, "--rdf_name_col_splitter", help="Divider when multiple annotations", type="character", default=",")
    parser <- add_argument(parser, "--db_name_col", help="Name for column containing ID in db", type="character")

    parser <- add_argument(parser, "--na_val", help="NA value field", default="NA", type="character")

    parser <- add_argument(parser, "--clean_up_descr_type", help="Clean up type, supported: TAIR", default=NA, type="character")
    parser <- add_argument(parser, "--clean_up_descr_colname", help="Description column", type="character", default=NA)
    
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
