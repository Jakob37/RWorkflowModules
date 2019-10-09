#!/usr/bin/env Rscript

library(argparser)

main <- function() {
    
    argv <- parse_input_params()
    
    suppressPackageStartupMessages(library(readr))

    allowed_types <- c("MaxQuant", "Proteios")
    if (!(argv$type %in% allowed_types)) {
        
    }
    
    if (argv$type == "Proteios") {
        rdf <- read_tsv(argv$in_fp, na="null", comment = "#")
        colnames(rdf) <- make.names()
    }
    else if (argv$type == "MaxQuant") {
        rdf <- read_tsv(argv$in_fp, na="0")
        colnames(rdf) <- gsub("^Intensity\\.", argv$prefix, make.names(colnames(rdf)))
    }
    else {
        stop("Provided type ", argv$type, " not in allowed types: ", paste(allowed_types, collapse=", "))
    }

    message("Writing output with dimensions: ", paste(dim(rdf), collapse=", "), " to: ", argv$out_fp)
    message("Header: ")
    print(colnames(rdf))
    write_tsv(rdf, path = argv$out_fp)    
}

parse_input_params <- function() {
    
    parser <- arg_parser("Quickly preparse LFQ expression matrices")
    parser <- add_argument(parser, "--in_fp", help="Input path", type="character", nargs=1)
    parser <- add_argument(parser, "--out_fp", help="Input path", type="character", nargs=1)
    parser <- add_argument(parser, "--type", help="Input path", type="character", nargs=1)
    parser <- add_argument(parser, "--prefix", help="Add prefix to sample names (MaxQuant)", type="character", default="")
    argv <- parse_args(parser)
    argv
}

if (!interactive()) {
    main()
}
