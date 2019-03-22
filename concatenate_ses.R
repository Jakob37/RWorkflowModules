#!/usr/bin/env Rscript

library(argparser)

parser <- arg_parser("Concatenate prepared ses objects into a nested and a flat master-object")
parser <- add_argument(parser, "--input_ses",    help="List of input ses objects to merge", type="character", nargs=Inf)
parser <- add_argument(parser, "--out_nested",   help="Path for nested (one list entry for each se)", type="character")
parser <- add_argument(parser, "--out_flat",     help="Path for flat (one list for all)", type="character")
parser <- add_argument(parser, "--trim_pattern", help="Trim this part of file to obtain sample pattern", type="character", default="_ses_obj.rds$")

argv <- parse_args(parser)

print(argv)

labels <- gsub(argv$trim_pattern, "", gsub(".*/", "", argv$input_ses))

print(labels)

flat_ses <- list()
nested_ses <- list()

for (i in seq_len(length(argv$input_ses))) {
    
    ses_path <- argv$input_ses[i]
    label <- labels[i]
    ses <- readRDS(ses_path)

    to_nested_ses <- ses
    nested_ses[[label]] <- to_nested_ses
    
    to_flat_ses <- ses
    names(to_flat_ses) <- paste(label, names(ses), sep="_")
    flat_ses <- c(flat_ses, to_flat_ses)
}

saveRDS(flat_ses, file=argv$out_flat)
saveRDS(nested_ses, file=argv$out_nested)
