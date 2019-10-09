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
    
    supported_cleanups <- c("TAIR", "SwissProt", "id_descr_os")
    
    if (argv$clean_up_descr_type %in% supported_cleanups) {
        if (is.na(argv$clean_up_descr_colname)) {
            message("No clean_up_colname assigned to clean up, skipping...")
        }
        else {
            message("Cleaning up annotations - might not be well adjusted for data! Trimming everything beyond inline-pipe character, and after 'OS='")
            if (argv$clean_up_descr_type == "TAIR") {
                db_df[[paste0(argv$clean_up_descr_colname, "_descr")]] <- db_df[[argv$clean_up_descr_colname]] %>% gsub("^\\| ", "", .) %>% gsub(" \\|.*", "", .) %>% gsub("OS=.*", "", .)
            }
            else if (argv$clean_up_descr_type == "SwissProt") {
                db_df[[paste0(argv$clean_up_descr_colname, "_descr")]] <- db_df[[argv$clean_up_descr_colname]] %>% gsub("^.* ", "", .)
            }
            else if (argv$clean_up_descr_type == "id_descr_os") {
                db_df[[paste0(argv$clean_up_descr_colname, "_descr")]] <- db_df[[argv$clean_up_descr_colname]] %>% 
                    gsub("\\w\\w\\|", "", .) %>%
                    gsub("\\|.*$", "", .)
            }
            else {
                stop("Unknown state, somethis is wrong in the code!")
            }
        }
    }
    else if (!is.na(argv$clean_up_descr_type)) {
        message("Only ", paste(supported_cleanups, collapse=", "), " supported as clean_up_annot type")
    }
    
    print(head(db_df[[argv$clean_up_descr_colname]]))
    
    protein_col <- argv$rdf_name_col
    
    rdf_annotated_protein_ids_list <- str_split(rdf[[protein_col]], argv$rdf_name_col_splitter)
    rdf_proteins_long <- tibble(
        id=seq_along(rdf_annotated_protein_ids_list), Protein=rdf_annotated_protein_ids_list
    ) %>% unnest() %>% data.frame()
    protein_annot_df <- rdf_proteins_long %>% dplyr::left_join(db_df, setNames(argv$db_name_col, "Protein"), keep=TRUE)

    combine_annotations <- function(df, group_info) {
        df %>% 
            filter(!is.na(Protein)) %>%
            dplyr::distinct(Protein, .keep_all = TRUE) %>%
            apply(2, function(col) { paste(col, collapse=",") }) %>% 
            t() %>% 
            data.frame() 
    }
    
    # Retrive unique matches per line, where multiple distinct are column-wise concatenated delimited by a comma    
    out <- protein_annot_df %>%
        dplyr::group_by(id) %>%
        dplyr::group_map(combine_annotations)
    
    # # Retrive unique matches per line, where multiple distinct are column-wise concatenated delimited by a comma    
    # out <- protein_annot_df %>%
    #     dplyr::group_by(id) %>%
    #     dplyr::group_map(
    #         function(df, group_info) { 
    #             df %>% 
    #                 filter(!is.na(UQ(as.name(protein_col)))) %>%
    #                 dplyr::distinct(UQ(as.name(protein_col)), .keep_all = TRUE) %>%
    #                 apply(2, function(col) { paste(col, collapse=",") }) %>% 
    #                 t() %>% 
    #                 data.frame() 
    #         })
    # )
    # 
    
    protein_annot_df <- do.call("rbind", out)
    protein_annot_df$ProteinIDFirst <- protein_annot_df[[protein_col]] %>% as.character() %>% strsplit(",") %>% lapply(function(entries) { entries[1] }) %>% unlist()
    
    message("Annotation count distribution - Number of entries with corresponding number matching proteins")
    print(protein_annot_df[[protein_col]] %>% unlist() %>% unname() %>% as.character() %>% strsplit(",") %>% lapply(., function(entries) { length(entries) }) %>% unlist() %>% table())
    
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
