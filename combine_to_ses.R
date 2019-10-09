#!/usr/bin/env Rscript

library(argparser)

parser <- arg_parser("Combine multiple datasets into a stored list of SummarizedExperiment objects")
parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character", nargs=Inf)
parser <- add_argument(parser, "--rdf_fps", help="Raw matrix path", type="character", nargs=Inf)
parser <- add_argument(parser, "--out_fp", help="Combined SummarizedExperiments file path", type="character")
parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character", nargs=Inf)
parser <- add_argument(parser, "--name", help="Name of the dataset in the list", default="")

argv <- parse_args(parser)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(SummarizedExperiment))

ddf <- readr::read_tsv(argv$ddf_fp, col_types = cols()) %>% data.frame()

rdfs <- lapply(argv$rdf_fps, function(rdf_fp) {
    readr::read_tsv(rdf_fp, col_types = cols()) %>% data.frame()
})

ses <- lapply(
    rdfs, 
    function(rdf, ddf, sample_col) {
        sdf <- rdf %>% select(one_of(ddf[[sample_col]])) %>% as.matrix()
        adf <- rdf %>% select(-one_of(ddf[[sample_col]]))
        SummarizedExperiment::SummarizedExperiment(
            assays=list(raw=sdf),
            rowData=adf,
            colData=ddf
        )
    },
    ddf=ddf,
    sample_col=argv$sample_col
)

if (argv$name == "") {
    names(ses) <- gsub("_[a-z]+_[a-z]+.tsv$", "", gsub(".*/", "", argv$rdf_fps))
} else {
    names(ses) <- argv$name
}
saveRDS(ses, file = argv$out_fp)
