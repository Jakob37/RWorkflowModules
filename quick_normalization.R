#!/usr/bin/env Rscript

library(argparser)

message("Requires NormalyzerDE to run (github.com/ComputationalProteomics/NormalyzerDE)")

parser <- arg_parser("Normalization wrapper for batch analysis")
parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
parser <- add_argument(parser, "--na_val", help="NA value fields", type="character", default="NA")

parser <- add_argument(parser, "--log2_fp", help="Output matrix path", type="character", default="")
parser <- add_argument(parser, "--median_fp", help="Output matrix path", type="character", default="")
parser <- add_argument(parser, "--loess_fp", help="Output matrix path", type="character", default="")
parser <- add_argument(parser, "--vsn_fp", help="Output matrix path", type="character", default="")

argv <- parse_args(parser)
# print(argv)

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

source("~/src/RWorkflowModules/debug_tools.R")
debug_tools$use_print_argv(argv)

message("Loaded raw data with dimensions: ", paste(dim(rdf), collapse=", "))

generate_normalized_data <- function(ddf, sdf, adf, norm_method, sample_col, out_fp) {
    norm_sdf <- norm_method(sdf)
    norm_rdf <- cbind(adf, norm_sdf)
    message("Writing ", nrow(norm_rdf), " rows to ", out_fp)
    readr::write_tsv(norm_rdf, path=out_fp)
}

if (argv$log2_fp != "") {
    generate_normalized_data(ddf, sdf, adf, base::log2, argv$sample_col, argv$log2_fp)
}

if (argv$median_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::medianNormalization, argv$sample_col, argv$median_fp)
}

if (argv$loess_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::performCyclicLoessNormalization, argv$sample_col, argv$loess_fp)
}

if (argv$vsn_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::performVSNNormalization, argv$sample_col, argv$vsn_fp)
}


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

generate_normalized_data <- function(ddf, sdf, adf, norm_method, sample_col, out_fp) {
    norm_sdf <- norm_method(sdf)
    norm_rdf <- cbind(adf, norm_sdf)
    message("Writing ", nrow(norm_rdf), " rows to ", out_fp)
    readr::write_tsv(norm_rdf, path=out_fp)
}

if (argv$log2_fp != "") {
    generate_normalized_data(ddf, sdf, adf, base::log2, argv$sample_col, argv$log2_fp)
}

if (argv$median_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::medianNormalization, argv$sample_col, argv$median_fp)
}

if (argv$loess_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::performCyclicLoessNormalization, argv$sample_col, argv$loess_fp)
}

if (argv$vsn_fp != "") {
    generate_normalized_data(ddf, sdf, adf, NormalyzerDE::performVSNNormalization, argv$sample_col, argv$vsn_fp)
}
