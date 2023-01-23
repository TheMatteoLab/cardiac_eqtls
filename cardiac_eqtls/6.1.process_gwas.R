setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/6.1.coloc_gwas"                             , showWarnings = FALSE)
dir.create("pipeline/6.1.coloc_gwas/coloc.gene"                  , showWarnings = FALSE)
dir.create("pipeline/6.1.coloc_gwas/coloc.isoform"               , showWarnings = FALSE)
#dir.create("pipeline/6.1.coloc_gwas/coloc.gene/coloc_by_trait"   , showWarnings = FALSE)
#dir.create("pipeline/6.1.coloc_gwas/coloc.isoform/coloc_by_trait", showWarnings = FALSE)

manifest                  = fread("input/gwas/manifest.txt", sep = "\t", header = TRUE, data.table = FALSE)
traits                    = readLines("input/gwas/pan_ukbb.filtered_traits.txt")
manifest                  = manifest[manifest$filename %in% traits, ]
manifest$n_controls_total = rowSums(manifest[, grepl("^n_controls", colnames(manifest)) == TRUE], na.rm = TRUE)
manifest$id               = gsub("\\.tsv.bgz$", "", manifest$filename)
manifest                  = manifest[, c("id", "trait_type", "phenocode", "description", "description_more", "coding_description", "category", "n_cases_full_cohort_both_sexes", "n_controls_total", "saige_heritability_EUR")]
manifest$filename         = paste(getwd(), "input/gwas/pan_ukbb/summary_statistics", paste(manifest$id, "txt", "gz", sep = "."), sep = "/")

fwrite(manifest, "pipeline/6.1.coloc_gwas/traits.manifest.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
