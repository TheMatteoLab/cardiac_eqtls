setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/4.1.differential_expression", showWarnings = FALSE)

# Input data

gene_info               = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE, data.table = FALSE)
isof_info               = fread("pipeline/1.2.expression//isoform_info.txt", sep = "\t", header = TRUE, data.table = FALSE)
gene_info$transcript_id = gene_info$gene_id

metadata   =              fread("pipeline/3.1.covariates/metadata.txt"     , sep = "\t", header = TRUE , data.table = FALSE)
covariates = add_rownames(fread("pipeline/3.1.covariates/covariates.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
int_list   = readLines         ("pipeline/3.2.eqtls/vars/cardiac_eqtls.txt")
cell_list  = readLines         ("pipeline/3.2.eqtls/vars/covariates_to_interaction.txt")
#int_list   = int_list[grepl("^peer", int_list) == FALSE & grepl("^pc", int_list) == FALSE] # original
int_list   = c("sex", "total_reads_norm") # remove mitochondrial reads as covariate
cell_list  = cell_list[grepl("^cibersort", cell_list) == TRUE]

covariates = covariates[,c(int_list, cell_list)]
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

run_diffexp = function(x, name, metadata, covariates, tissue1, tissue2, int_list, cell_list, gene_info)
{
	message(paste(name, tissue1, tissue2))
	
	runs              = intersect(metadata$run, rownames(covariates))
	metadata          = metadata  [metadata$run %in% runs,]
	covariates        = covariates[runs,]
	covariates        = covariates[metadata[metadata$tissue %in% c(tissue1, tissue2), "run"],]
	covariates$tissue = 0
	
	covariates[metadata[metadata$tissue == tissue1, "run"], "tissue"] = 1
	
	x                 = x[, rownames(covariates)]
	diffexp           = as.data.frame(rbindlist(lapply(rownames(x), function(gene_id)
	{
		out = as.data.frame(rbindlist(lapply(cell_list, function(cell)
		{
			this_covariates = covariates[,c("tissue", int_list, cell)]
			this_out        = run_diffexp_lm(gene_id, x, this_covariates)
			this_out$cell   = cell
			
			return(this_out)
		})), stringsAsFactors = FALSE)
		
		out$transcript_id = gene_id
		
		return(out)
	})), stringsAsFactors = FALSE)
	
	
	#rownames(diffexp) = rownames(x)
	#diffexp$qval      = qvalue(diffexp$pval)$qvalues
	diffexp$qval      = p.adjust(diffexp$pval, method = "bonferroni") 
	diffexp$tissue1   = tissue1
	diffexp$tissue2   = tissue2
	diffexp           = merge(gene_info[,c("transcript_id", "gene_id", "gene_name", "gene_type")], diffexp)
	
	fwrite(diffexp, paste("pipeline/4.1.differential_expression", paste("diffexp_tissue_by_cell", name, tissue1, tissue2, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

invisible(lapply(c("ipsc_cvpc", "heart"), function(tissue1)
{
	lapply(c("heart", "arteria"), function(tissue2)
	{
		if(tissue1 != tissue2)
		{
			#run_diffexp(gene_tpm, "gene_tpm"   , metadata, covariates, tissue1, tissue2)
			#run_diffexp(isof_tpm, "isoform_tpm", metadata, covariates, tissue1, tissue2)
			run_diffexp(isof_use, "isoform_use", metadata, covariates, tissue1, tissue2, int_list, cell_list, isof_info)
		}
	})
}))



