setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/4.1.differential_expression", showWarnings = FALSE)

# Input data
metadata   =              fread("pipeline/3.1.covariates/metadata.txt"     , sep = "\t", header = TRUE , data.table = FALSE)
covariates = add_rownames(fread("pipeline/3.1.covariates/covariates.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
int_list   = readLines         ("pipeline/3.2.eqtls/vars/cardiac_eqtls.gene.txt")
int_cell   = readLines         ("pipeline/3.2.eqtls/vars/covariates_to_interaction.txt")
#int_list   = int_list[grepl("^peer", int_list) == FALSE & grepl("^pc", int_list) == FALSE] # original
int_list   = c("sex", "total_reads_norm") # remove mitochondrial reads as covariate
covariates = covariates[,c(int_list, int_cell)]
gene_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_gene.normalized.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
isof_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))
isof_use   = add_rownames(fread("pipeline/1.2.expression/use_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))

# Run differential expression vs. CIBERSORT cell types

run_diffexp_lm = function(gene_id, x, covariates, vars0, var1)
{
	covariates$exp    = as.numeric(x[gene_id, rownames(covariates)])
	mylm              = as.data.frame(coefficients(summary(lm(exp ~ ., data = covariates[, c("exp", var1, vars0)]))))
	out               = mylm[var1, ]
	colnames(out)     = c("beta", "se", "tval", "pval")
	out$cell_type     = var1
	out$transcript_id = gene_id
	
	return(out)
}

run_diffexp = function(x, name, metadata, covariates, tissue1, tissue2)
{
	message(name)
	
	runs              = intersect(metadata$run, rownames(covariates))
	metadata          = metadata  [metadata$run %in% runs,]
	covariates        = covariates[runs,]
	x                 = x[,rownames(covariates)]
	vars0             = c("sex", "total_reads_norm", "uniquely_mapped_reads_to_canonical_chromsomes", "mitochondrial_reads", 'ipsc_cvpc', 'heart', 'arteria') # original
	vars0             = c("sex", "total_reads_norm", 'ipsc_cvpc', 'heart', 'arteria')# remove mitochondrial reads as covariate
	vars1             = colnames(covariates)[grepl("cibersort" , colnames(covariates)) == TRUE]
	diffexp           = as.data.frame(rbindlist(lapply(rownames(x), function(gene_id)
	{
		as.data.frame(rbindlist(lapply(vars1, function(var1)
		{
			run_diffexp_lm(gene_id, x, covariates, vars0, var1)
		})), stringsAsFactors = FALSE)
	})), stringsAsFactors = FALSE)
	
	diffexp$qval      = p.adjust(diffexp$pval, method = "bonferroni") 
	
	fwrite(diffexp, paste("pipeline/4.1.differential_expression", paste("diffexp_cell", name, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

run_diffexp(gene_tpm, "gene_tpm"   , metadata, covariates)
run_diffexp(isof_tpm, "isoform_tpm", metadata, covariates)
run_diffexp(isof_use, "isoform_use", metadata, covariates)
