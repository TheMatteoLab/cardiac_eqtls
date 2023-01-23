
setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

library(coloc)

dir.create(paste("pipeline/3.2.eqtls", "eqtl_overlap", "cardiac_eqtls.coloc_by_gene"   , sep = "/"), showWarnings = FALSE)

geneinfo_genes    = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE , data.table = FALSE)
geneinfo_isoforms = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
eqtl_genes        = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.gene.egenes.txt"     , sep = "\t", header = TRUE, data.table = FALSE)
eqtl_isoforms     = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.isoform.egenes.txt"  , sep = "\t", header = TRUE, data.table = FALSE)
eqtl_isoforms     = eqtl_isoforms[eqtl_isoforms$egene == TRUE,]

intersect_overlapping_genes = function(gene_id, name, eqtls, eqtl_genes)
{
    this       = unique(eqtls[eqtls$gene_id == gene_id, "transcript_id"])
    eqtls      = eqtls[eqtls$transcript_id %in% this,]
    this       = c(this, gene_id)
    eqtl_genes = eqtl_genes[eqtl_genes$gene_id == gene_id, ]
    
    if(nrow(eqtl_genes) > 1){eqtl_genes = eqtl_genes[eqtl_genes$egene == TRUE,]}
    
    eqtls    = rbind(eqtls, eqtl_genes)
    analysis = "isoform"
    coloc    = as.data.frame(rbindlist(lapply(1:(length(this) - 1), function(ii)
    {
        transcript_id1 = this[[ii]]
        eqtl1          = fread(paste("pipeline/3.2.eqtls", "eqtls_by_gene", paste(name, analysis, sep = "."), paste("qtl", transcript_id1, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
        types1         = sort(unique(eqtls[eqtls$transcript_id == transcript_id1, "type"]))
        
        out = as.data.frame(rbindlist(lapply(types1, function(type1)
        {
            eqtl_type1 = eqtl1[eqtl1$type == type1,]
            out = as.data.frame(rbindlist(lapply((ii + 1):length(this), function(jj)
            {
                analysis2      = "isoform"
                transcript_id2 = this[[jj]]
                if(grepl("ENSG", transcript_id2) == TRUE){analysis2 = "gene"}
				
				infile2        = paste("pipeline/3.2.eqtls", "eqtls_by_gene", paste(name, analysis2, sep = "."), paste("qtl", transcript_id2, "txt", sep = "."), sep = "/")
				
				if(file.exists(infile2) == TRUE)
				{
					eqtl2          = fread(infile2, sep = "\t", header = TRUE, data.table = FALSE)
					types2         = sort(unique(eqtls[eqtls$transcript_id == transcript_id2, "type"]))

					out = as.data.frame(rbindlist(lapply(types2, function(type2)
					{
						eqtl_type2 = eqtl2[eqtl2$type == type2,]
						tocoloc    = merge(eqtl_type1[,c("id", "pval", "af")], eqtl_type2[,c("id", "pval")], by = "id", suffixes = 1:2)

						if(nrow(tocoloc) > 100)
						{
							coloc_mapped = coloc.abf(dataset1 = list(snp = tocoloc$id, pvalues = tocoloc$pval1, N = 955, MAF = tocoloc$af, type = "quant"),
													 dataset2 = list(snp = tocoloc$id, pvalues = tocoloc$pval2, N = 955, MAF = tocoloc$af, type = "quant")
													) 

							probs           = as.data.frame(t(coloc_mapped$summary))
							myres           = coloc_mapped$results
							myres           = myres[, c(which(colnames(myres) == "snp"), ncol(myres))]
							colnames(myres) = c("id", "pp_snp")
							myres           = merge(tocoloc, myres)
							myres           = merge(unique(eqtl_type1[,c("id", "chrom", "pos", "ref", "alt", "af")]), myres)
							myres           = cbind(data.frame(transcript_id1 = transcript_id1, 
															   transcript_id2 = transcript_id2, 
															   type1          = type1, 
															   type2          = type2) , myres)
							myres           = myres[order(myres$pp_snp, decreasing = TRUE), ]
							out             = cbind(probs, myres[1, ])

							fwrite(myres, paste("pipeline/3.2.eqtls/eqtl_overlap", 
												paste(name, "coloc_by_gene", sep = "."), 
												paste(transcript_id1, type1, transcript_id2, type2, "txt", sep = "."), 
												sep = "/"), 
								   sep = "\t", col.names = TRUE, row.names = FALSE)

						}else
						{
							out = data.frame(nsnps          = nrow(tocoloc),
											 PP.H0.abf      = 1,
											 PP.H1.abf      = 0,
											 PP.H2.abf      = 0,
											 PP.H3.abf      = 0,
											 PP.H4.abf      = 0,
											 transcript_id1 = transcript_id1, 
											 transcript_id2 = transcript_id2, 
											 type1          = type1, 
											 type2          = type2
											)
							out = cbind(out, eqtl_type1[which.min(eqtl_type1$pval), c("id", "af", "chrom", "pos", "ref", "alt")])
							out = cbind(out, data.frame(pval1 = min(eqtl_type1$pval), pval2 = 1, pp_snp = 0))
						}
						return(out)
					})), stringsAsFactors = FALSE)
				}else
				{
					out = data.frame(nsnps          = nrow(eqtl_type1),
									 PP.H0.abf      = 1,
									 PP.H1.abf      = 0,
									 PP.H2.abf      = 0,
									 PP.H3.abf      = 0,
									 PP.H4.abf      = 0,
									 transcript_id1 = transcript_id1, 
									 transcript_id2 = transcript_id2, 
									 type1          = type1, 
									 type2          = 0
									)
					out = cbind(out, eqtl_type1[which.min(eqtl_type1$pval), c("id", "af", "chrom", "pos", "ref", "alt")])
					out = cbind(out, data.frame(pval1 = min(eqtl_type1$pval), pval2 = 1, pp_snp = 0))
				}
            })), stringsAsFactors = FALSE)
        })), stringsAsFactors = FALSE)
    })), stringsAsFactors = FALSE)
    
    return(coloc)
}

intersect_qtls = function(name, eqtls, eqtl_genes)
{
    transcript_ids = sort(unique(eqtls$transcript_id))
    gene_ids       = sort(unique(eqtls$gene_id      ))
    
    out     = as.data.frame(rbindlist(lapply(gene_ids, function(gene_id){intersect_overlapping_genes(gene_id, name, eqtls, eqtl_genes)})), stringsAsFactors = FALSE)
    
    fwrite(out, paste("pipeline/3.2.eqtls/eqtl_overlap", 
                      paste(name, "coloc_by_gene", "txt", sep = "."), 
                      sep = "/"), 
           sep = "\t", col.names = TRUE, row.names = FALSE)
    
    return(out)
}

coloc = intersect_qtls("cardiac_eqtls", eqtl_isoforms, eqtl_genes)
