setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

option_list   = list(make_option("--name"            , type="character", default="test1", help="output folder name", metavar="character"),
					 make_option("--analysis"        , type="character", default="gene" , help="gene or isoform"   , metavar="character")
					) 
					
opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
name              = opt$name
analysis          = opt$analysis

#name              = "cardiac_eqtls"
#analysis          = "gene"

if(analysis == "gene"   ){geneinfo = fread("pipeline/1.2.expression/gene_info.txt"   , sep = "\t", header = TRUE , data.table = FALSE)}
if(analysis == "isoform"){geneinfo = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)}

folder = paste(name, analysis, sep = ".")
qtls   = merge_qtls(folder, geneinfo)

message("####################################################################################")
message(paste("Tested genes"      , length(unique(qtls$transcript_id)), sep = " = "))
message(paste("eGenes"            , paste(nrow(qtls[qtls$egene == TRUE & qtls$type == 0, ]), " (", signif(nrow(qtls[qtls$egene == TRUE & qtls$type == 0, ]) / nrow(qtls[                     qtls$type == 0, ]) * 100, digits = 3), "%)", sep = ""), sep = " = "))
message(paste("Conditional eGenes", paste(nrow(qtls[qtls$egene == TRUE & qtls$type == 1, ]), " (", signif(nrow(qtls[qtls$egene == TRUE & qtls$type == 1, ]) / nrow(qtls[qtls$egene == TRUE & qtls$type == 0, ]) * 100, digits = 3), "%)", sep = ""), sep = " = "))


merge_ints = function(infolder, qtls)
{
    infiles                = list.files(paste("pipeline/3.2.eqtls//eqtls_by_gene", infolder, sep = "/"), pattern = "^int", full.names = TRUE)
    indata                 = suppressWarnings(as.data.frame(rbindlist(lapply(infiles, function(x){fread(x, sep = "\t", header = TRUE, data.table = FALSE)})), stringsAsFactors = FALSE))
	colnames(indata)       = c("chrom", "pos", "ref", "alt", "rsid", "id", "af", "beta", "se", "pval", "transcript_id", "bonferroni", "fdr", "tests", "type", "beta_int", "se_int", "pval_int", "interaction", "beta_0", "se_0", "pval_0", "beta_1", "se_1", "pval_1")
	indata$transcript2type = paste(indata$transcript_id, indata$type)
	qtls  $transcript2type = paste(qtls  $transcript_id, qtls  $type)
	transcript2types       = qtls[qtls$egene == TRUE, "transcript2type"]
	indata                 = indata[indata$transcript2type %in% transcript2types,]
	indata$transcript2type = NULL
	out                    = merge(unique(qtls[,c("transcript_id", "gene_id", "gene_name", "gene_type", "start", "end", "strand")]), indata)
    out$distance           = out$pos - out$start

    out[out$strand == "-", "distance"] = out[out$strand == "-", "end"] - out[out$strand == "-", "pos"]
	
	interactions = sort(unique(out$interaction))
	interactions = interactions[grepl("cibersort.combined", interactions) == FALSE & grepl("ipsc_cvpc.cibersort", interactions) == FALSE & grepl("adult.cibersort", interactions) == FALSE]
	
	fdr = as.data.frame(rbindlist(lapply(interactions, function(x)
	{
		this            = out[out$interaction == x,]
		this$qval_int   = p.adjust(this$pval_int, method = "bonferroni")
		this$qval_1     = p.adjust(this$pval_1  , method = "BH")
		this$qval_0     = p.adjust(this$pval_0  , method = "BH")
		this$int_signif = FALSE

		this[this$qval_int <= 0.05, "int_signif"] = TRUE
		
		this$cell       = FALSE
		this$specific   = FALSE
		this$associated = FALSE
		
		if(x == "ipsc_cvpc"){this_adult = this}

		this[this$int_signif == TRUE & this$qval_1 < 0.05 & this$qval_0 < 0.05 & abs(this$beta_1) > abs(this$beta_0), "associated"] = TRUE
		this[this$int_signif == TRUE & this$qval_1 < 0.05 & this$qval_0 > 0.05                                      , "specific"  ] = TRUE
		
		this[this$associated == TRUE | this$specific == TRUE, "cell"] = TRUE
		
		if(x == "ipsc_cvpc")
		{
			this_adult$interaction = "adult"
			this_adult[, c("beta_0", "se_0", "pval_0", "qval_0", "beta_1", "se_1", "pval_1", "qval_1")] = this_adult[, c("beta_1", "se_1", "pval_1", "qval_1", "beta_0", "se_0", "pval_0", "qval_0")]
			this_adult[this_adult$int_signif == TRUE & this_adult$qval_1 < 0.05 & this_adult$qval_0 < 0.05 & abs(this_adult$beta_1) > abs(this_adult$beta_0), "associated"] = TRUE
			this_adult[this_adult$int_signif == TRUE & this_adult$qval_1 < 0.05 & this_adult$qval_0 > 0.05                                                  , "specific"  ] = TRUE
			
			this_adult[this_adult$associated == TRUE | this_adult$specific == TRUE, "cell"] = TRUE
			
			this = rbind(this, this_adult)
		}
		
		return(this)
	})), stringsAsFactors = FALSE)
    
    fwrite(fdr, paste("pipeline/3.2.eqtls/eqtls", paste(infolder, "interactions", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
    
    return(fdr)
}

ints = merge_ints(folder, qtls)

message("")
message(paste("Genes with tissue or cell type associations", paste(length(unique(ints[ints$cell == TRUE, "transcript_id"])), " (", signif(length(unique(ints[ints$cell == TRUE, "transcript_id"])) / nrow(qtls[qtls$egene == TRUE & qtls$type == 0, ]) * 100, digits = 3), "%)", sep = ""), sep = " = "))

invisible(lapply(sort(unique(ints[grepl("cibersort", ints$interaction) == TRUE, "interaction"])), function(x)
{
	this = ints[ints$interaction == x & ints$cell == TRUE,]
	
	message(paste(gsub("_", " ", sub("^cibersort.regular\\.", "", x)), " = ", length(unique(this$transcript_id)), " (Specific = ", length(unique(this[this$specific == TRUE, "transcript_id"])), "; Associated = ", length(unique(this[this$associated == TRUE, "transcript_id"])), ")", sep = ""))
}))

invisible(lapply(sort(unique(ints[grepl("cibersort", ints$interaction) == FALSE, "interaction"])), function(x)
{
	this = ints[ints$interaction == x & ints$cell == TRUE,]
	
	message(paste(gsub("_", "-", x), " = ", length(unique(this$transcript_id)), " (Specific = ", length(unique(this[this$specific == TRUE, "transcript_id"])), "; Associated = ", length(unique(this[this$associated == TRUE, "transcript_id"])), ")", sep = ""))
}))

message("####################################################################################")
message("")

