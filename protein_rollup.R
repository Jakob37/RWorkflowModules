#!/usr/bin/env Rscript

library(argparser)

parser <- arg_parser("Statistical comparison wrapper")
parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
parser <- add_argument(parser, "--out_fp", help="Output matrix path", type="character")

parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
parser <- add_argument(parser, "--peptide_col", help="Peptide column in main data frame", type="character")
parser <- add_argument(parser, "--protein_col", help="Protein column in main data frame", type="character")

parser <- add_argument(parser, "--protein_tools_path", help="CraftOmics protein tools path", type="character")

argv <- parse_args(parser)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(limma))

source(argv$protein_tools_path)

pr$protein_rollup_on_matrix(
    design_fp=argv$ddf_fp,
    data_fp=argv$rdf_fp,
    protein_col=argv$protein_col,
    peptide_col=argv$peptide_col,
    sample_col=argv$sample_col,
    out_fp=argv$out_fp,
    rollup_func="rrollup"
)
