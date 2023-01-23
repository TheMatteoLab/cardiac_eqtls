### Version history:
### V02: 6/10/2020
### V01: 4/17/2020

# General functions

add_rownames = function(x) # add rownames to fread
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

### Gene expression normalization
my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
    df                           = as.matrix(df)
    data_valid_expressed_full_qn = normalize.quantiles(df, copy=FALSE)
    input_mat                    = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}

normalize_tpm = function(outfolder, tpm, geneinfo, normalize = TRUE, gene_ids = NULL, filter_exp = FALSE, min_tpm = 2, min_samples = 0.1, use = FALSE)
{
    if(is.null(gene_ids) == FALSE){tpm = tpm[gene_ids,]}
    message("Normalizing expression...")
    message(paste("Total genes/peaks", nrow(tpm), sep = " = "))
    message(paste("Total samples"    , ncol(tpm), sep = " = "))
    
    expressed = as.matrix(tpm)
    
    if(filter_exp == TRUE)
    {
        expressed[as.matrix(tpm) <  min_tpm] = 0
        expressed[as.matrix(tpm) >= min_tpm] = 1

        tpm_f = tpm[rowSums(expressed) >= (min_samples * ncol(tpm)),]
        
    }else
    {
        tpm_f = tpm
    }
    
    if(use == TRUE)
    {
        geneinfo   = geneinfo[geneinfo$transcript_id %in% rownames(tpm_f),]
		gene_table = table(geneinfo$gene_id)
		geneinfo   = geneinfo[geneinfo$gene_id %in% names(gene_table[gene_table > 1]),]
		
        tpm_f = as.data.frame(rbindlist(lapply(sort(unique(geneinfo$gene_id)), function(gene_id)
        {
            this_transcript_ids = geneinfo[geneinfo$gene_id == gene_id, "transcript_id"]
            this                = tpm_f[this_transcript_ids, ]
            mysums              = unlist(lapply(colSums(this), function(x){max(c(100, x))}))
            out                 = as.data.frame(100 * t(t(as.matrix(this)) / mysums))
            out$transcript_id   = this_transcript_ids
            
            return(out)
        })), stringsAsFactors = FALSE)
        
        rownames(tpm_f)     = tpm_f$transcript_id
        tpm_f$transcript_id = NULL
    }
    
    message(paste("Total expressed genes"    , nrow(tpm_f), sep = " = "))
    
    if(normalize == TRUE)
    {
        tpm_f_std_norm = transform_standard_normal(tpm_f)
    }else
    {
        tpm_f_std_norm = tpm_f
    }
    
    fwrite    (tpm_f          , file = paste(outfolder, "expressed" , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite    (tpm_f_std_norm , file = paste(outfolder, "normalized", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    writeLines(rownames(tpm_f), con  = paste(outfolder, "gene_ids"  , "txt", sep = "."))
    
    return(list(tpm_f = tpm_f, tpm_f_std_norm = tpm_f_std_norm, gene_ids = rownames(tpm_f), sample_ids = colnames(tpm_f)))
}

calculate_peer_factors = function(outfolder, expdata, peer_factor_n = 10, n_genes = 0)
{
    tpm_f_std_norm    = expdata$tpm_f_std_norm
    tpm_f             = expdata$tpm_f
    
    if (n_genes > 0)
    {
        totest         = data.frame(gene_id = rownames(tpm_f), sd = unlist(apply(tpm_f, 1, sd)))
        totest         = totest[order(totest$sd, decreasing = TRUE), "gene_id"]
        tpm_f_std_norm = tpm_f_std_norm[totest[1:n_genes],]
    }
   
    model = PEER()
    PEER_setPhenoMean (model, t(as.matrix(tpm_f_std_norm)))
    PEER_setNk        (model, peer_factor_n               )
    PEER_update       (model                              )

    factors   = PEER_getX        (model)
    weights   = PEER_getW        (model)
    precision = PEER_getAlpha    (model)
    residuals = PEER_getResiduals(model)
    precision = data.frame(peer = paste("peer", 1:peer_factor_n, sep = "") , precision = precision, stringsAsFactors = FALSE)
    residuals = t(residuals)

    rownames(factors  ) = colnames(tpm_f_std_norm)
    rownames(weights  ) = rownames(tpm_f_std_norm)
    rownames(residuals) = rownames(tpm_f_std_norm)
    colnames(residuals) = colnames(tpm_f_std_norm)
    colnames(factors  ) = paste("peer", 1:peer_factor_n, sep = "")
    colnames(weights  ) = paste("peer", 1:peer_factor_n, sep = "")

    factors          = as.data.frame(factors  )
    weights          = as.data.frame(weights  )
    residuals        = as.data.frame(residuals)
    factors$assay_id = rownames     (factors  )
    
    fwrite(factors  , file = paste(outfolder, "peer", "factors"  , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(weights  , file = paste(outfolder, "peer", "weights"  , "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(precision, file = paste(outfolder, "peer", "precision", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    fwrite(residuals, file = paste(outfolder, "peer", "residuals", "txt", sep = "."), sep = "\t", row.names = TRUE , col.names = TRUE)
    
    return(factors)
}

divide_phenotypes_by_gene = function(expdata, outfolder)
{
    message("Dividing phenotype data by gene/peak...")
    
    dir.create(outfolder, showWarnings = FALSE)

    gene_ids   = expdata$gene_ids
    sample_ids = expdata$sample_ids
    rawexp     = expdata$tpm_f
    normexp    = expdata$tpm_f_std_norm
    
    invisible(lapply(gene_ids, function(gene)
    {
        outdata = data.frame(sample_id = sample_ids, raw = as.numeric(rawexp[gene,sample_ids]), norm = as.numeric(normexp[gene,sample_ids]))
        fwrite(outdata, paste(outfolder, paste(gene, "txt", sep = "."), sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
    }))
}

divide_phenotypes_tissue = function(outfolder, tpm, geneinfo, run_peer = FALSE, normalize = FALSE, gene_ids = NULL, filter_exp = FALSE, phenotype_min_value = 1, phenotype_min_samples = 0.1, peer_factor_n = 10, n_genes = 1000, use = FALSE)
{
    message("---------------------------")
    
    expdata = normalize_tpm(outfolder, tpm, geneinfo, normalize, gene_ids, filter_exp, min_tpm = phenotype_min_value, min_samples = phenotype_min_samples, use)
    
    if(run_peer == TRUE){peerdata = calculate_peer_factors(outfolder, expdata, peer_factor_n, n_genes)}
    
    divide_phenotypes_by_gene(expdata, outfolder)
}

# eQTL analysis
write_h5_file = function(h5_file, analysis, gtinfo, gtdata, gene_id, expdata)
{
    invisible(suppressWarnings(file.remove(h5_file)))

    invisible(h5createFile (h5_file))
    invisible(h5createGroup(h5_file, "genotype"))
    invisible(h5createGroup(h5_file, "genotype/col_header"))
    invisible(h5createGroup(h5_file, "genotype/row_header"))
    invisible(h5createGroup(h5_file, "phenotype"))
    invisible(h5createGroup(h5_file, "phenotype/col_header"))
    invisible(h5createGroup(h5_file, "phenotype/row_header"))
    
    genedata = geneinfo[geneinfo$transcript_id == gene_id,]

    h5write(gtinfo$chrom           , file = h5_file, name="genotype/col_header/chrom"  )
    h5write(gtinfo$pos             , file = h5_file, name="genotype/col_header/pos"    )
    h5write(gtinfo$pos             , file = h5_file, name="genotype/col_header/pos_cum")
    h5write(as.matrix(gtdata)      , file = h5_file, name="genotype/matrix"            )
    h5write(colnames (gtdata)      , file = h5_file, name="genotype/row_header/sample_ID")
    h5write(c(gene_id)             , file = h5_file, name="phenotype/col_header/gene_ID"     )
    h5write(genedata[, "chrom"    ], file = h5_file, name="phenotype/col_header/gene_chrom"  )
    h5write(genedata[, "end"      ], file = h5_file, name="phenotype/col_header/gene_end"    )
    h5write(genedata[, "start"    ], file = h5_file, name="phenotype/col_header/gene_start"  )
    h5write(genedata[, "strand"   ], file = h5_file, name="phenotype/col_header/gene_strand" )
    h5write(genedata[, "gene_name"], file = h5_file, name="phenotype/col_header/phenotype_ID")
    h5write(as.matrix(expdata)     , file = h5_file, name="phenotype/matrix"                 )
    h5write(colnames( expdata)     , file = h5_file, name="phenotype/row_header/sample_ID")
}

run_eigenmt = function(gt_file, geneloc_file, snploc_file, qtl_file, fdr_file, chrom)
{
    command = paste("python", paste(getwd(),"script/eigenMT.py", sep = "/"), 
                    "--CHROM"   , chrom,
                     "--QTL"    , qtl_file,
                     "--GEN"    , gt_file,
                     "--GENPOS" , snploc_file,
                     "--PHEPOS" , geneloc_file,
                     "--OUT"    , fdr_file
                   )
				   
    system(command)
}


run_eqtl = function(gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
{
    suppressWarnings(write_h5_file(h5_file, analysis, gtinfo, gtdata, gene_id, expdata))
	
    command_py = paste("python", paste(getwd(), "script", "run_limix.py", sep = "/"), tmp_folder, gene_id, "scan", "normal")
    
    system(command_py)
    
    indata                   = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
	colnames(indata)         = c("beta", "se", "pval")
    indata                   = cbind(gtinfo, indata)
    indata$gene_id           = gene_id
    indata$bonferroni        = p.adjust(indata$pval, method = "bonferroni")
	
	if(nrow(indata[indata$se > 100, ]) > 0){indata[indata$se > 100, "se"] = 100}
	
    if(min(indata$bonferroni) == 1)
    {
        fdrdata       = indata[which.min(indata$pval), ]
        fdrdata$fdr   = 1
        fdrdata$tests = nrow(indata)
    }else
    {
        eigenmt_input            = indata[,c("id", "gene_id", "beta", "se", "pval", "bonferroni")]
        colnames(eigenmt_input)  = c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
        qtl_file_emt_in          = sub("qtl.csv", "eigenmt_input.txt" , qtl_file_tmp)
        qtl_file_emt_out         = sub("qtl.csv", "eigenmt_output.txt", qtl_file_tmp)
	
        fwrite(eigenmt_input, qtl_file_emt_in , sep = "\t", row.names = FALSE, col.names = TRUE)
	
        run_eigenmt(gt_file, geneloc_file, snploc_file, qtl_file_emt_in, qtl_file_emt_out, chrom)
	
        fdrdata           = fread(qtl_file_emt_out, sep = "\t", header = TRUE, data.table = FALSE)
        colnames(fdrdata) = c("id", "gene_id", "beta", "se", "pval", "bonferroni", "fdr", "tests")
    }
    
    fdrdata           = merge(gtinfo, fdrdata)
    fdrdata           = fdrdata[,c(colnames(indata), "fdr", "tests")]
    return(list(lead = fdrdata, qtl = indata))
	return(indata)
}



# Find interactions
find_interactions = function(totest, var1, gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
{
    suppressWarnings(write_h5_file(h5_file, analysis, gtinfo, gtdata, gene_id, expdata))
	
    command_py = paste("python", paste(getwd(), "script", "run_limix.py", sep = "/"), tmp_folder, gene_id, var1, "normal")
    
    system(command_py)
    
    indata             = add_rownames(fread(qtl_file_tmp, sep = ",", header = TRUE, data.table = FALSE))
	colnames(indata)   = c("beta", "se", "pval")
	indata             = cbind(totest, indata)
	indata$interaction = var1
	
	return(indata)
}

# Merge eQTLs

merge_qtls = function(infolder, geneinfo)
{
    infiles = list.files(paste("pipeline/3.2.eqtls//eqtls_by_gene", infolder, sep = "/"), pattern = "^fdr", full.names = TRUE)
    
    indata     = suppressWarnings(as.data.frame(rbindlist(lapply(infiles, function(x){fread(x, sep = "\t", header = TRUE, data.table = FALSE)})), stringsAsFactors = FALSE))
    outdata    = list()
    gene_ids   = unique(indata$gene_id)
    
    for(type in sort(unique(indata$type)))
    {
        if(length(gene_ids) > 0)
        {
            this = indata[indata$type == type & indata$gene_id %in% gene_ids,]

            #this$qval  = qvalue(this$fdr)$qvalues
            this$qval  = p.adjust(this$fdr, method = "BH")
            this$egene = FALSE

            this[this$qval <= 0.05, "egene"] = TRUE

            gene_ids            = this[this$egene == TRUE, "gene_id"]
            outdata[[type + 1]] = this
        }else
        {
            break
        }
    }
    
    out          = as.data.frame(rbindlist(outdata), stringsAsFactors = FALSE)
    out          = merge(geneinfo[,c("transcript_id", "gene_id", "gene_name", "gene_type", "start", "end", "strand")], out, by.x = "transcript_id", by.y = "gene_id")
    out$distance = out$pos - out$start

    out[out$strand == "-", "distance"] = out[out$strand == "-", "end"] - out[out$strand == "-", "pos"]
    
    fwrite(out, paste("pipeline/3.2.eqtls/eqtls", paste(infolder, "egenes", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
    
    return(out)
}


