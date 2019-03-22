#!/usr/bin/env Rscript

library(argparser)

parser <- arg_parser("Statistical comparison wrapper")
parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
parser <- add_argument(parser, "--out_fp", help="Output matrix path", type="character")
parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character")
parser <- add_argument(parser, "--na_val", help="NA value fields", type="character", default=NA)
parser <- add_argument(parser, "--omit_samples", help="List of samples to omit", type="character", nargs=Inf)

parser <- add_argument(parser, "--model", help="Statistical model to use", type="character")
parser <- add_argument(parser, "--contrasts", help="Contrasts to perform", type="character", nargs = Inf)

argv <- parse_args(parser)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(limma))

# source("~/src/RWorkflowModules/debug_tools.R")
# debug_tools$use_print_argv(argv)

ddf <- readr::read_tsv(argv$ddf_fp, col_types=cols())
rdf <- readr::read_tsv(argv$rdf_fp, col_types=cols(), comment = "#", na = argv$na_val)
sdf <- rdf %>% dplyr::select(ddf[[argv$sample_col]]) %>% as.matrix()
adf <- rdf %>% dplyr::select(-one_of(ddf[[argv$sample_col]]))

if (!any(is.na(argv$omit_samples))) {
    stat_sdf <- sdf %>% as.data.frame() %>% dplyr::select(-one_of(argv$omit_samples)) %>% as.matrix()
} else {
    stat_sdf <- sdf
}

message("Loaded raw data with dimensions: ", paste(dim(rdf), collapse=", "))

# Prepare R equation objects
# factor_ddf <- data.frame(do.call("cbind", lapply(data.frame(ddf), function(col) { as.factor(as.character(col)) })))

model <- as.formula(argv$model)
model_design <- model.matrix(model, data=ddf)

# Calculate Limma table (pairwise)
fit <- limma::lmFit(sdf, model_design)
contrast.matrix <- limma::makeContrasts(contrasts=argv$contrasts, levels=model_design)

fit_contrasts <- contrasts.fit(fit, contrast.matrix)
fit_bayes <- limma::eBayes(fit_contrasts)

limma_tables <- lapply(
    seq_len(length(colnames(fit_bayes$coefficients))), 
    function(coef) { 
        my_tbl <- topTable(fit_bayes, coef=coef, number=Inf, sort.by="none")
        my_tbl$absLogFC <- abs(my_tbl$logFC)
        my_tbl
    })
names(limma_tables) <- colnames(fit_bayes$coefficients)
combined_limma_tables <- do.call("cbind", limma_tables)
combined_limma_tables <- cbind(row_nbr=seq_len(nrow(combined_limma_tables)), combined_limma_tables)

# Calculate Limma table (F-statistic), draft
# fit_bayes_fstats <- limma::eBayes(fit)
# f_stat_table <- limma::topTable(fit_bayes_fstats, number=Inf, sort.by="none")
# colnames(f_stat_table) <- paste("fstat", colnames(f_stat_table), sep="_")

# Output
full_stat_rdf <- cbind(adf, combined_limma_tables, sdf)

message("Writing data with dimensions: ", paste(dim(full_stat_rdf), collapse=", "))
readr::write_tsv(full_stat_rdf, path=argv$out_fp)
