setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/4.15.differential_expression_ridge", showWarnings = FALSE)

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(hdi))

# Input data
metadata   =              fread("pipeline/3.1.covariates/metadata.txt"     , sep = "\t", header = TRUE , data.table = FALSE)
covariates = add_rownames(fread("pipeline/3.1.covariates/covariates.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
int_list   = readLines         ("pipeline/3.2.eqtls/vars/cardiac_eqtls.txt")
int_list   = c(int_list[grepl("^peer", int_list) == FALSE & grepl("^pc", int_list) == FALSE], colnames(covariates)[grepl("cibersort.regular", colnames(covariates)) == TRUE])
covariates = covariates[,int_list]
gene_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_gene.normalized.txt"   , sep = "\t", header = TRUE , data.table = FALSE))
isof_tpm   = add_rownames(fread("pipeline/1.2.expression/tpm_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))
isof_use   = add_rownames(fread("pipeline/1.2.expression/use_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))

# Analysis
run_diffexp_ridge_lm = function(gene_id, x, covariates)
{
	covariates$exp = as.numeric(x[gene_id, rownames(covariates)])
    
    mylm = ridge.proj(x = as.matrix(covariates[, colnames(covariates) != "exp"]), y = covariates$exp, suppress.grouptesting = FALSE)
    out  = data.frame(transcript_id = gene_id, covariate = names(mylm$pval), beta = mylm$bhat, se = mylm$se, pval = mylm$pval)
    
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
	
	x                 = x[, rownames(covariates)]
	diffexp           = as.data.frame(rbindlist(lapply(rownames(x), function(gene_id){run_diffexp_ridge_lm(gene_id, x, covariates)})), stringsAsFactors = FALSE)
	diffexp$qval      = p.adjust(diffexp$pval, method = "bonferroni") 
	diffexp$tissue1   = tissue1
	diffexp$tissue2   = tissue2
	
	fwrite(diffexp, paste("pipeline/4.15.differential_expression_ridge", paste("diffexp", name, tissue1, tissue2, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

#invisible(lapply(c("ipsc_cvpc", "heart"), function(tissue1)
#{
#	lapply(c("heart", "arteria"), function(tissue2)
#	{
#		if(tissue1 != tissue2)
#		{
#			run_diffexp(gene_tpm, "gene_tpm"   , metadata, covariates, tissue1, tissue2)
#			#run_diffexp(isof_tpm, "isoform_tpm", metadata, covariates, tissue1, tissue2)
#			run_diffexp(isof_use, "isoform_use", metadata, covariates, tissue1, tissue2)
#		}
#	})
#}))
#
gene_info               = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE, data.table = FALSE)
isof_info               = fread("pipeline/1.2.expression//isoform_info.txt", sep = "\t", header = TRUE, data.table = FALSE)
gene_info$transcript_id = gene_info$gene_id

read_diffexp = function(name, tissue1, tissue2, gene_info)
{
    indata               = fread(paste("pipeline/4.15.differential_expression_ridge", paste("diffexp", name, tissue1, tissue2, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
    #indata$transcript_id = rownames(indata)
    indata               = indata[order(indata$pval),]
    indata               = merge(gene_info[,c("transcript_id", "gene_id", "gene_name", "gene_type")], indata)
    indata$qval          = p.adjust(indata$pval / length(unique(indata$covariate)), method = "bonferroni")
    indata$tissue1       = tissue1
    indata$tissue2       = tissue2
    indata$type          = name
    indata$diffexp       = FALSE
    
    indata[indata$qval < 0.05, "diffexp"] = TRUE
    
    message(paste(name, tissue1, tissue2, nrow(indata), nrow(indata[indata$qval < 0.05,])))
    return(indata)
}

diffexp = as.data.frame(rbindlist(lapply(c("ipsc_cvpc", "heart"), function(tissue1)
{
	out = as.data.frame(rbindlist(lapply(c("heart", "arteria"), function(tissue2)
	{
		if(tissue1 != tissue2)
		{
			indata_gene    = read_diffexp("gene_tpm"   , tissue1, tissue2, gene_info)
			#indata_iso_tpm = read_diffexp("isoform_tpm", tissue1, tissue2, isof_info)
			indata_iso_use = read_diffexp("isoform_use", tissue1, tissue2, isof_info)
            
            #return(rbind(indata_gene, indata_iso_tpm, indata_iso_use))
            return(rbind(indata_gene, indata_iso_use))
		}
	})), stringsAsFactors = FALSE)
    
     return(out)
})), stringsAsFactors = FALSE)

#name    = "gene_tpm"
#tissue1 = "ipsc_cvpc"
#tissue2 = "heart"
#head(read_diffexp("gene_tpm"   , tissue1, tissue2, gene_info))
#head(read_diffexp("isoform_use", tissue1, tissue2))

fwrite(diffexp, "pipeline/4.15.differential_expression_ridge/diffexp.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

