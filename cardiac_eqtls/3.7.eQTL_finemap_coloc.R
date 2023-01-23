
setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

library(coloc)

dir.create(paste("pipeline/3.2.eqtls", "eqtls_fine_map", "cardiac_eqtls.gene"   , sep = "/"), showWarnings = FALSE)
dir.create(paste("pipeline/3.2.eqtls", "eqtls_fine_map", "cardiac_eqtls.isoform", sep = "/"), showWarnings = FALSE)

geneinfo_genes    = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE , data.table = FALSE)
geneinfo_isoforms = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
eqtl_genes        = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.gene.egenes.txt"         , sep = "\t", header = TRUE, data.table = FALSE)
eqtl_isoforms     = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.isoform.egenes.txt"      , sep = "\t", header = TRUE, data.table = FALSE)
eqtl_genes        = eqtl_genes   [eqtl_genes   $egene == TRUE,]
eqtl_isoforms     = eqtl_isoforms[eqtl_isoforms$egene == TRUE,]
covariates        = add_rownames(fread    ("pipeline/3.1.covariates/covariates.txt", sep = "\t", header = TRUE , data.table = FALSE))
n                 = nrow(covariates)

quiet = function(x) 
{ 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 

finemap_qtl = function(name, analysis, transcript_id, type, n, geneinfo)
{
    eqtl              = fread(paste("pipeline/3.2.eqtls", "eqtls_by_gene", paste(name, analysis, sep = "."), paste("qtl", transcript_id, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
    this              = eqtl[eqtl$type == type,]
    finemap           = finemap.abf(dataset=list(pvalues = this$pval, MAF = this$af, snp = this$id, N = n, sdY = 1, type = "quant"))
    finemap           = finemap[,c("snp", "pvalues.", "SNP.PP")]
    colnames(finemap) = c("id", "pval", "pp")
    
	finemap             = finemap[finemap$id != "null",]
	finemap$pp          = finemap$pp/sum(finemap$pp)
	finemap             = finemap[order(finemap$pp, -finemap$pval, decreasing = TRUE),]
    finemap$cum         = cumsum(finemap$pp)
    finemap$cs          = FALSE
    to_cs               = finemap[finemap$cum >= 0.99, ][1, "id"]
    
    finemap[finemap$cum < 0.99 | finemap$id == to_cs, "cs"] = TRUE

    fwrite(finemap, paste("pipeline/3.2.eqtls", "eqtls_fine_map", paste(name, analysis, sep = "."), paste("qtl", transcript_id, type, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)

    credible_set = finemap[ finemap$cs == TRUE, "id"]
	
	out               = finemap[finemap$id %in% credible_set,]
	out$cs            = NULL
	out$transcript_id = transcript_id
	out$type          = type
    
    #out = merge(eqtl[,c("transcript_id", "type")], finemap[finemap$id %in% credible_set,])
    #out = merge(geneinfo[,c("transcript_id", "gene_id", "gene_name")], out, by.x = "transcript_id", by.y = "gene_id")
    
    return(out[order(out$pp, decreasing = TRUE),])
}

finemap_genes    = as.data.frame(rbindlist(lapply(1:nrow(eqtl_genes   ), function(ii){quiet(finemap_qtl("cardiac_eqtls", "gene"   , eqtl_genes   [ii, "transcript_id"], eqtl_genes   [ii, "type"], n, geneinfo_genes   ))})), stringsAsFactors = FALSE)
finemap_isoforms = as.data.frame(rbindlist(lapply(1:nrow(eqtl_isoforms), function(ii){quiet(finemap_qtl("cardiac_eqtls", "isoform", eqtl_isoforms[ii, "transcript_id"], eqtl_isoforms[ii, "type"], n, geneinfo_isoforms))})), stringsAsFactors = FALSE)

fwrite(finemap_genes   , paste("pipeline/3.2.eqtls", "eqtls_fine_map", paste("cardiac_eqtls", "gene"   , "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(finemap_isoforms, paste("pipeline/3.2.eqtls", "eqtls_fine_map", paste("cardiac_eqtls", "isoform", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
