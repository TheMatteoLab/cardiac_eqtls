setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

option_list   = list(make_option("--taskid"          , type="integer"  , default=0      , help="SGE task ID"       , metavar="character")) 

opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
taskid            = opt$taskid

exp_all             = add_rownames(fread("pipeline/1.2.expression//use_isoform.normalized.txt", sep = "\t", header = TRUE , data.table = FALSE))
gt_folder           = paste(getwd(), "pipeline", "1.3.genotype"  , "tpm_gene"   , sep = "/")
geneinfo_gene       = fread("pipeline/1.2.expression/gene_info.txt"                          , sep = "\t", header = TRUE , data.table = FALSE)
geneinfo_isoform    = fread("pipeline/1.2.expression/isoform_info.txt"                       , sep = "\t", header = TRUE , data.table = FALSE)
eqtl_rbp            = fread("pipeline/7.1.trans/rbp/input/eqtl_rbp.txt"                      , sep = "\t", header = TRUE , data.table = FALSE)
eqtl                = eqtl_rbp[taskid, c("gene_id", "gene_name", "type", "chrom", "pos", "ref", "alt", "id", "beta", "se", "pval", "fdr", "qval", "egene")]
rbp                 = eqtl$gene_id
totest              = fread(paste("pipeline/7.1.trans/rbp/input/rbp2transcript", paste(rbp, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
vars0               = readLines(paste("pipeline/3.2.eqtls", "vars", paste("cardiac_eqtls", "isoform", "txt", sep = "."), sep = "/"))
exp_rbp             = add_rownames(fread(paste(paste(getwd(), "pipeline", "1.2.expression", "tpm_gene", sep = "/"), paste(rbp, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE))
covariates          = add_rownames(fread("pipeline/3.1.covariates/covariates.txt", sep = "\t", header = TRUE , data.table = FALSE))
metadata            =              fread("pipeline/3.1.covariates/metadata.txt"  , sep = "\t", header = TRUE , data.table = FALSE)
gtdata              = add_rownames(fread(paste(gt_folder , paste("gt_data", rbp, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE))
covariates          =              merge(metadata, covariates, by.x = "run", by.y = "row.names")
covariates          = merge(covariates, data.frame(wgs_id = colnames(gtdata), gt = as.numeric(gtdata[eqtl$id,])))

calculate_trans = function(x, name, tolm, vars0)
{
    tolm$totest    = tolm[,x]
    out0           = as.data.frame(summary(lm(norm_isoform ~ ., data = tolm[,c("norm_isoform", "totest", "norm", vars0)]))$coefficients)
    out0           = out0["totest", c("Estimate", "Std. Error", "Pr(>|t|)")]
    colnames(out0) = paste(c("beta", "se", "pval"), name, sep = "_")
    
    return(out0)
}

out = as.data.frame(rbindlist(lapply(rownames(exp_all), function(isoform)
{
    expdata                  = data.frame(run = colnames(exp_all), norm_isoform = as.numeric(exp_all[isoform,]))
    tolm                     = merge(covariates, expdata, by   = "run")
    tolm                     = merge(tolm      , exp_rbp, by.x = "run", by.y = "row.names")
    out_gt                   = calculate_trans("gt"  , "gt" , tolm, vars0)
    out_ex                   = calculate_trans("norm", "exp", tolm, vars0)
    out0                     = cbind(eqtl[,c("gene_id", "gene_name", "type")], out_gt, out_ex)
    out0$transcript_id_trans = isoform  
    
    return(out0)
})), stringsAsFactors = FALSE)

fwrite(out, paste("pipeline/7.1.trans/rbp/trans_eqtls_by_rbp", paste(rbp, eqtl$type, "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)

