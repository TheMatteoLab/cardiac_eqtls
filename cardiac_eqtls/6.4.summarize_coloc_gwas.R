setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

option_list   = list(make_option("--analysis"        , type="character", default="gene" , help="gene or isoform"   , metavar="character")) 
					
opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
analysis          = opt$analysis

if(analysis == "gene"   ){geneinfo = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE , data.table = FALSE)}
if(analysis == "isoform"){geneinfo = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)}

eqtls = fread(paste("pipeline/3.2.eqtls/eqtls/cardiac_eqtls", analysis, "egenes.txt"      , sep = "."), sep = "\t", header = TRUE , data.table = FALSE)
ints  = fread(paste("pipeline/3.2.eqtls/eqtls/cardiac_eqtls", analysis, "interactions.txt", sep = "."), sep = "\t", header = TRUE , data.table = FALSE)

infolder             = paste("pipeline/6.1.coloc_gwas/coloc", analysis, sep = ".")
egenes               = unique(eqtls[eqtls$egene == TRUE, c("gene_id", "gene_name", "transcript_id")])
egenes$infile        = paste(infolder, paste("top", egenes$transcript_id, "txt", sep = "."), sep = "/")
egenes$infile_exists = unlist(lapply(egenes$infile, file.exists))

message("####################################################################################")
message(paste("Tested genes"      , length(unique(egenes$gene_id      )), sep = " = "))
message(paste("Tested isoforms"   , length(unique(egenes$transcript_id)), sep = " = "))

coloc = as.data.frame(rbindlist(lapply(egenes$infile, function(x)
{
	indata = fread(x, sep = "\t", header = TRUE , data.table = FALSE)
	return(indata)
})), stringsAsFactors = FALSE)

coloc     = coloc[coloc$nsnps > 0,]
coloc$tr2type = paste(coloc$transcript_id, coloc$type)
ints $tr2type = paste(ints $transcript_id, ints $type)

coloc2int = merge(coloc, ints[ints$cell == TRUE,], by = "tr2type", suffixes = c("_coloc", "_eqtl"))

fwrite(coloc    , paste("pipeline/6.1.coloc_gwas/coloc.eqtls"       , analysis, "txt", sep = "."), sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(coloc2int, paste("pipeline/6.1.coloc_gwas/coloc.interactions", analysis, "txt", sep = "."), sep = "\t", col.names = TRUE, row.names = FALSE)

message(paste("eQTLs that colocalize with GWAS"                     , length(unique(coloc    [coloc    $PP.H4.abf > 0.5, "tr2type"])), sep = " = "))
message(paste("Stage-associated eQTLs that colocalize with GWAS"    , length(unique(coloc2int[coloc2int$PP.H4.abf > 0.5 & coloc2int$interaction %in% c("ipsc_cvpc", "adult"), "tr2type"])), sep = " = "))
message(paste("Tissue-associated eQTLs that colocalize with GWAS"   , length(unique(coloc2int[coloc2int$PP.H4.abf > 0.5 & coloc2int$interaction %in% c("heart", "arteria", "heart_atrium", "heart_ventricle", "arteria_aorta", "arteria_coronary"), "tr2type"])), sep = " = "))
message(paste("Cell type-associated eQTLs that colocalize with GWAS", length(unique(coloc2int[coloc2int$PP.H4.abf > 0.5 & grepl("^cibersort", coloc2int$interaction) == TRUE, "tr2type"])), sep = " = "))

