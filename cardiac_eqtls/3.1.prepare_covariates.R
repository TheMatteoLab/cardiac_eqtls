setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/3.1.covariates", showWarnings = FALSE)

covariates_rna             = fread("pipeline/1.2.expression/covariates_rna.txt"  , sep = "\t", header = TRUE , data.table = FALSE)
metadata                   = fread("pipeline/1.2.expression/metadata.txt"        , sep = "\t", header = TRUE , data.table = FALSE)
meta_subject               = fread("pipeline/1.2.expression/meta_subject.txt"    , sep = "\t", header = TRUE , data.table = FALSE)
meta_rna                   = fread("pipeline/1.1.metadata/meta_rna.txt"          , sep = "\t", header = TRUE , data.table = FALSE)
pca_data                   = fread("input/genotypes/genotype_pca.txt"            , sep = " " , header = FALSE, data.table = FALSE)
cibersort_output           = fread("pipeline/2.1.scrna_seq/cibersort_results.txt", sep = "\t", header = TRUE , data.table = FALSE)
pca_data                   = add_rownames(add_rownames(pca_data        ))
cibersort_output           =              add_rownames(cibersort_output)
colnames(pca_data)         = paste("pc", 1:ncol(pca_data), sep = "")
colnames(cibersort_output) = paste("cibersort", colnames(cibersort_output), sep = ".")

colnames(covariates_rna)[[1]] = "run"

covariates = metadata
covariates = merge(covariates, pca_data        , by.x = "wgs_id", by.y = "row.names")
covariates = merge(covariates, covariates_rna  )
covariates = merge(covariates, cibersort_output, by.x = "run"   , by.y = "row.names")
covariates = merge(covariates, meta_subject    )
covariates = merge(covariates, meta_rna[,c("run", "tissue", "body_site")])

covariates$total_reads_norm = covariates$total_reads / mean(covariates$total_reads)


rownames(covariates) = covariates$run

meta2qtl = covariates[,c("run", "assay_id", "wgs_id", "subject_id", "subject_name", "study", "tissue", "body_site")]
covs2qtl = covariates[,c("sex", "age", "height", "weight", "bmi", "total_reads", "total_reads_norm", "uniquely_mapped_reads_to_canonical_chromsomes", "mitochondrial_reads", paste("peer", 1:300, sep = ""), paste("pc", 1:20, sep = ""), colnames(cibersort_output))]

for(x in sort(unique(meta2qtl$body_site)))
{
	covs2qtl[,x] = 0
	covs2qtl[meta2qtl[meta2qtl$body_site == x, "run"], x] = 1
}

for(x in c("heart", "arteria"))
{
	covs2qtl[,x] = 0
	covs2qtl[meta2qtl[meta2qtl$tissue == x, "run"], x] = 1
}

fwrite(meta2qtl, "pipeline/3.1.covariates/metadata.txt"  , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(covs2qtl, "pipeline/3.1.covariates/covariates.txt", sep = "\t", col.names = TRUE, row.names = TRUE )

# rewrite kinship

kinship     = fread("pipeline/1.4.kinship/kinship.rel"   , sep = "\t", header = FALSE, data.table = FALSE)
kinship_ids = fread("pipeline/1.4.kinship/kinship.rel.id", sep = "\t", header = FALSE, data.table = FALSE)[,1]

rownames(kinship) = kinship_ids
colnames(kinship) = kinship_ids

fwrite(kinship, "pipeline/3.1.covariates/kinship.txt", sep = "\t", col.names = TRUE, row.names = TRUE )

# calculate average % cells/tissue

cell2tissue    = as.data.frame(t(100 * as.matrix(add_rownames(aggregate(. ~ tissue   , data = covariates[, grepl("cibersort.regular", colnames(covariates)) == TRUE | colnames(covariates) %in% c("tissue"   )], FUN = mean)))))
cell2body_site = as.data.frame(t(100 * as.matrix(add_rownames(aggregate(. ~ body_site, data = covariates[, grepl("cibersort.regular", colnames(covariates)) == TRUE | colnames(covariates) %in% c("body_site")], FUN = mean)))))

cell2body_site$ipsc_cvpc = NULL

cell2tissue = cbind(cell2tissue, cell2body_site)

fwrite(cell2tissue, "pipeline/3.1.covariates/cell2tissue.txt", sep = "\t", col.names = TRUE, row.names = TRUE )
