setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/4.1.differential_expression", showWarnings = FALSE)

# Input data
metadata   =              fread("pipeline/3.1.covariates/metadata.txt"     , sep = "\t", header = TRUE , data.table = FALSE)
covariates = add_rownames(fread("pipeline/3.1.covariates/covariates.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
int_list   = readLines         ("pipeline/3.2.eqtls/vars/cardiac_eqtls.txt")
#int_list   = int_list[grepl("^peer", int_list) == FALSE & grepl("^pc", int_list) == FALSE] # original
int_list   = c("sex", "total_reads_norm") # remove mitochondrial reads as covariate
covariates = covariates[,int_list]
gene_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_gene.normalized.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
isof_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))
isof_use   = add_rownames(fread("pipeline/1.2.expression/use_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))

# Analysis
run_diffexp_lm = function(gene_id, x, covariates)
{
	covariates$exp = as.numeric(x[gene_id, rownames(covariates)])
	mylm           = as.data.frame(coefficients(summary(lm(exp ~ ., data = covariates))))
	out            = mylm["tissue", ]
	colnames(out)  = c("beta", "se", "tval", "pval")
	
	return(out)
}

run_diffexp = function(x, name, metadata, covariates, tissue1, tissue2)
{
	message(paste(name, tissue1, tissue2))
	
	runs              = intersect(metadata$run, rownames(covariates))
	metadata          = metadata  [metadata$run %in% runs,]
	covariates        = covariates[runs,]
	covariates        = covariates[metadata[metadata$tissue %in% c(tissue1, tissue2), "run"],]
	covariates$tissue = 0
	
	covariates[metadata[metadata$tissue == tissue1, "run"], "tissue"] = 1
	
	x                 = x[,rownames(covariates)]
	diffexp           = as.data.frame(rbindlist(lapply(rownames(x), function(gene_id){run_diffexp_lm(gene_id, x, covariates)})), stringsAsFactors = FALSE)
	rownames(diffexp) = rownames(x)
	#diffexp$qval      = qvalue(diffexp$pval)$qvalues
	diffexp$qval      = p.adjust(diffexp$pval, method = "bonferroni") 
	diffexp$tissue1   = tissue1
	diffexp$tissue2   = tissue2
	
	fwrite(diffexp, paste("pipeline/4.1.differential_expression", paste("diffexp", name, tissue1, tissue2, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = TRUE)
}

invisible(lapply(c("ipsc_cvpc", "heart"), function(tissue1)
{
	lapply(c("heart", "arteria"), function(tissue2)
	{
		if(tissue1 != tissue2)
		{
			run_diffexp(gene_tpm, "gene_tpm"   , metadata, covariates, tissue1, tissue2)
			run_diffexp(isof_tpm, "isoform_tpm", metadata, covariates, tissue1, tissue2)
			run_diffexp(isof_use, "isoform_use", metadata, covariates, tissue1, tissue2)
		}
	})
}))



