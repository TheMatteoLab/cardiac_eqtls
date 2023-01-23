setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

suppressPackageStartupMessages(library(coloc))

option_list   = list(make_option("--taskid"          , type="integer"  , default=0      , help="SGE task ID"       , metavar="character"),
					 make_option("--analysis"        , type="character", default="gene" , help="gene or isoform"   , metavar="character")
					) 
					
opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
taskid            = opt$taskid
analysis          = opt$analysis

if(analysis == "gene")
{
	geneinfo   = fread("pipeline/1.2.expression/gene_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
	exp_folder = paste(getwd(), "pipeline", "1.2.expression", "tpm_gene", sep = "/")
	gt_folder  = paste(getwd(), "pipeline", "1.3.genotype"  , "tpm_gene", sep = "/")
	out_folder = paste(getwd(), "pipeline", "3.2.eqtls"     , "tpm_gene", sep = "/")
}
if(analysis == "isoform")
{
	geneinfo   = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
	exp_folder = paste(getwd(), "pipeline", "1.2.expression", "use_isoform", sep = "/")
	gt_folder  = paste(getwd(), "pipeline", "1.3.genotype"  , "use_isoform", sep = "/")
	out_folder = paste(getwd(), "pipeline", "3.2.eqtls"     , "use_isoform", sep = "/")
}

rownames(geneinfo) = geneinfo$transcript_id

chromsizes  = read.table("/frazer01/reference/public/hg19/hg19.size.txt", sep = "\t", header = FALSE, col.names = c("chrom", "size"))[,1:2]
gene_id     = geneinfo[taskid, "transcript_id"]
#gene_id     = "ENSG00000124588.20_5"
geneinfo    = geneinfo[geneinfo$transcript_id == gene_id,]
manifest    = fread("pipeline/6.1.coloc_gwas/traits.manifest.txt", sep = "\t", header = TRUE, data.table = FALSE)
manifest    = manifest[file.exists(manifest$filename) == TRUE,]
qtl_folder  = paste(getwd(), "pipeline/3.2.eqtls"     , "eqtls_by_gene", paste("cardiac_eqtls", analysis, sep = "."), sep = "/")
out_folder  = paste(getwd(), "pipeline/6.1.coloc_gwas",                  paste("coloc"        , analysis, sep = "."), sep = "/")
eqtls       = fread(paste("pipeline", "3.2.eqtls", "eqtls", paste("cardiac_eqtls", analysis, "egenes.txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
this_qtl    = eqtls[eqtls$transcript_id == gene_id & eqtls$egene == TRUE,]

#manifest = manifest[1:3,]

if(nrow(this_qtl) > 0)
{
	sink_outputs  = file(paste(out_folder, paste("TMP", gene_id, "sink", "output" , "txt", sep = "."), sep = "/"), open = "wt")
	sink(sink_outputs , type = "output")
	
	qtl_file   = paste(qtl_folder, paste("qtl", gene_id      , "txt", sep = "."), sep = "/")
	qtl_data   = fread(qtl_file, sep = "\t", header = TRUE, data.table = FALSE)
	coord      = paste(sub("chr", "", geneinfo$chrom), ":", max(c(0, geneinfo$start - 500000)), "-", geneinfo$end + 500000, sep = "")

	gwas_out    = lapply(manifest$id, function(gwas)
	{
		tmp_file   = paste(out_folder, paste("TMP", gene_id, gwas, "txt", sep = "."), sep = "/")
		gwas_file  = manifest[manifest$id == gwas, "filename"  ]
		trait_type = manifest[manifest$id == gwas, "trait_type"]
		n          = rowSums(manifest[manifest$id == gwas, c("n_cases_full_cohort_both_sexes", "n_controls_total")])
		command    = paste("tabix", gwas_file, coord, ">", tmp_file)

		system(command)
		
		if(file.size(tmp_file) > 0)
		{

			indata           = fread(tmp_file, sep = "\t", header = FALSE , data.table = FALSE)
			myhead           = unlist(strsplit(system(paste("zcat", gwas_file, "|", "head", "-n", 1), intern = TRUE), split = "\t"))
			colnames(indata) = myhead
			indata  $idx     = paste(indata  $chr  , indata  $pos, sep = "_")
			qtl_data$idx     = paste(qtl_data$chrom, qtl_data$pos, sep = "_")
			totest           = merge(qtl_data, indata, by = "idx", suffixes = 1:2)
			totest           = totest[ totest$pval_heterogeneity > 1e-6,]
			
			if(nrow(totest) > 0)
			{
				type2coloc       = "quant"

				if(!trait_type %in% c("biomarkers", "continuous"))
				{
					totest$af_meta = totest$af_controls_meta
					type2coloc     = "cc"
					cases_fr       = manifest[manifest$id == gwas, "n_cases_full_cohort_both_sexes"] / n
				}

				if(nrow(totest[totest$ref1 == totest$alt2 & totest$alt1 == totest$ref2,]) > 0)
				{
					totest1 = totest[totest$ref1 == totest$ref2 & totest$alt1 == totest$alt2,]
					totest2 = totest[totest$ref1 == totest$alt2 & totest$alt1 == totest$ref2,]
					
					totest2$beta_meta =   - totest2$beta_meta
					totest2$af_meta   = 1 - totest2$af_meta
					
					totest = rbind(totest1, totest2)
				}

				totest  = totest[order(totest$pos1),]
				outdata = lapply(this_qtl$type, function(type)
				{
					tocoloc = totest[totest$type == type,]
					
					if(trait_type %in% c("biomarkers", "continuous"))
					{
						coloc_mapped = coloc.abf(dataset1 = list(snp = tocoloc$id, pvalues = tocoloc$pval     , N = 966, MAF = tocoloc$af     , type = "quant"),
												 dataset2 = list(snp = tocoloc$id, pvalues = tocoloc$pval_meta, N = n  , MAF = tocoloc$af_meta, type = "quant"))
					}else
					{
						coloc_mapped = coloc.abf(dataset1 = list(snp = tocoloc$id, pvalues = tocoloc$pval     , N = 966, MAF = tocoloc$af     , type = "quant"),
												 dataset2 = list(snp = tocoloc$id, pvalues = tocoloc$pval_meta, N = n  , MAF = tocoloc$af_meta, type = "cc"   , s = cases_fr))
					}

					probs           = as.data.frame(t(coloc_mapped$summary))
					myres           = coloc_mapped$results
					myres           = myres[, c(which(colnames(myres) == "snp"), ncol(myres))]
					colnames(myres) = c("id", "pp_snp")
					myres           = merge(unique(qtl_data[,c("id", "chrom", "pos", "ref", "alt", "af")]), myres)
					myres           = cbind(data.frame(transcript_id = gene_id, type = type), myres)
					myres           = myres[order(myres$pp_snp, decreasing = TRUE), ]
					out             = cbind(probs, myres[1, ])
					out$trait       = gwas
					myres$trait     = gwas
					
					return(list(top = out, snp = myres))
				})
                
                outdata = list(top = as.data.frame(rbindlist(lapply(1:length(outdata), function(ii){outdata[[ii]][["top"]]})), stringsAsFactors = FALSE),
                               snp = as.data.frame(rbindlist(lapply(1:length(outdata), function(ii){outdata[[ii]][["snp"]]})), stringsAsFactors = FALSE)
                              )
                return(outdata)
			}else
			{
				out = data.frame(nsnps          = 0,
								 PP.H0.abf      = 1,
								 PP.H1.abf      = 0,
								 PP.H2.abf      = 0,
								 PP.H3.abf      = 0,
								 PP.H4.abf      = 0,
								 transcript_id  = gene_id, 
								 type           = this_qtl$type
								)
				out        = cbind(out, qtl_data[which.min(qtl_data$pval), c("id", "chrom", "pos", "ref", "alt", "af")])
				out$pp_snp = 0
				out$trait  = gwas
				myres      = data.frame(transcript_id  = gene_id, 
									    type           = this_qtl$type,
										id             = "",
										chrom          = "chr0",
										pos            = 0,
										ref            = "",
										alt            = "",
										af             = 0,
										pp_snp         = 0,
										trait          = gwas)

				outdata = list(top = out, snp = myres)
			}
		}else
		{
			out = data.frame(nsnps          = 0,
							 PP.H0.abf      = 1,
							 PP.H1.abf      = 0,
							 PP.H2.abf      = 0,
							 PP.H3.abf      = 0,
							 PP.H4.abf      = 0,
							 transcript_id  = gene_id, 
							 type           = this_qtl$type
							)
			out        = cbind(out, qtl_data[which.min(qtl_data$pval), c("id", "chrom", "pos", "ref", "alt", "af")])
			out$pp_snp = 0
			out$trait       = gwas
			
			myres      = data.frame(transcript_id  = gene_id, 
									type           = this_qtl$type,
									id             = "",
									chrom          = "chr0",
									pos            = 0,
									ref            = "",
									alt            = "",
									af             = 0,
									pp_snp         = 0,
									trait          = gwas)

			outdata = list(top = out, snp = myres)
		}
		
		return(outdata)
	})
	
	towrite_top = suppressWarnings(as.data.frame(rbindlist(lapply(1:length(gwas_out), function(ii){gwas_out[[ii]][["top"]]})), stringsAsFactors = FALSE))
	towrite_snp = suppressWarnings(as.data.frame(rbindlist(lapply(1:length(gwas_out), function(ii){gwas_out[[ii]][["snp"]]})), stringsAsFactors = FALSE))
	
	fwrite(towrite_top, paste(out_folder, paste("top", gene_id, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
	fwrite(towrite_snp, paste(out_folder, paste("snp", gene_id, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
	
	sink()
	
	system(paste("rm", paste(out_folder, paste("TMP", gene_id, "*", sep = "."), sep = "/")))
	#message(paste("rm", paste(out_folder, paste("TMP", gene_id, "*", sep = "."), sep = "/")))
}
