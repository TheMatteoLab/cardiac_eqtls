setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

option_list = list(make_option("--taskid", type="integer", default=0, help="SGE task ID", metavar="character")) 

opt_parser    = OptionParser(option_list=option_list)
opt           = parse_args(opt_parser)
taskid        = opt$taskid
qtl_distance  = 500000
outfolder     = paste(getwd(), "pipeline/1.3.genotype/tpm_gene", sep = "/")
maf_threshold = 0.05

filter_vcf_by_region = function(wgs_list, region, maf_threshold, vcf_input, outfile)
{
    if(grepl("ipscore", vcf_input) == TRUE)
    {
        region  = paste("chr", region, sep = "")
        tmpfile = sub("\\.gz$", "", outfile)
        outdata = paste("-O v", "-o", tmpfile)
    }else
    {
        outdata = paste("-O z", "-o", outfile)
    }
    command1 = paste("bcftools", "view", 
                     "-f", "PASS", 
                     "-r", region, 
                     "-q", paste(maf_threshold, "minor", sep = ":"), 
                     "-s", paste(wgs_list, collapse = ","), 
					 "--force-samples",
                     vcf_input, 
                     outdata
                    )
    
    command2 = paste("bcftools", "index", "--tbi", "-f", outfile)
    
    system(command1)
    if(grepl("ipscore", vcf_input) == TRUE)
    {
        command3 = paste("sed", '"s/^chr//"', tmpfile, "|", "bcftools", "view", "-O z", "-o", outfile)
        
        system(command3)
    }
    
    system(command2)
    
    return(outfile)
}

read_gttable = function(infile, wgs_list)
{
    gttable                = fread(infile, header = TRUE, sep = "\t", data.table = FALSE)
    colnames(gttable)      =  gsub(":GT", "", unlist(lapply(colnames(gttable), function(x){unlist(strsplit(x, "]"))[[2]]})))
	colnames(gttable)[1:6] = c("chrom", "pos", "ref", "alt", "rsid", "af")
    gttable$id        = paste("VAR", gttable$chrom, gttable$pos, gttable$ref, gttable$alt, sep = "_") 
    gttable$rsid      = gttable$rsid
    rownames(gttable) = gttable$id
    return(gttable)
}

chromsizes  = read.table("/frazer01/reference/public/hg19/hg19.size.txt", sep = "\t", header = FALSE, col.names = c("chrom", "size"))[,1:2]
metadata    = fread     ("pipeline/1.2.expression/metadata.txt"         , sep = "\t", header = TRUE , data.table = FALSE)
geneinfo    = fread     ("pipeline/1.2.expression/gene_info.txt"        , sep = "\t", header = TRUE , data.table = FALSE)
gene_id     = geneinfo[taskid, "gene_id"]
tmp_folder  = paste("/scratch", gene_id, sep = "/")
wgs_list    = sort(unique(metadata[metadata$vcf == TRUE, "wgs_id"]))
chrom       = sub("chr", "", geneinfo[geneinfo$gene_id == gene_id, "chrom"])
gene_from   =                geneinfo[geneinfo$gene_id == gene_id, "start"]
gene_to     =                geneinfo[geneinfo$gene_id == gene_id, "end"  ]
region_from = max(0                                            , gene_from - qtl_distance)
region_to   = min(chromsizes[chromsizes$chrom == chrom, "size"], gene_to   + qtl_distance)
region      = paste(chrom, ":", region_from, "-", region_to, sep = "")
study2files = data.frame(study = c("ipscore", "gtex"), vcf = c(vcf_file_ipscore, vcf_file_gtex))

dir.create(tmp_folder, showWarnings = FALSE)

vcf_list    = unlist(lapply(1:nrow(study2files), function(ii)
{
	this_wgs_list = sort(unique(metadata[metadata$study  == study2files[ii, "study"], "wgs_id"]))
	tmp_vcf       = paste(tmp_folder, paste("tmp", gene_id, study2files[ii, "study"], "vcf", "gz", sep = "."), sep = "/")
	tmp_vcf       = filter_vcf_by_region(this_wgs_list, region, maf_threshold, study2files[ii, "vcf"], tmp_vcf)
	return(tmp_vcf)
}))

tmp_vcf     = paste(tmp_folder, paste("tmp", gene_id, "all", "vcf", "gz", sep = "."), sep = "/")
command1    = paste("bcftools", "merge" , "-m", "none", paste(rev(vcf_list), collapse = " "), "|", 
					"bcftools", "norm"  , "-m-", "|", 
					"bcftools", "filter", "-i", "'F_PASS(GT!=\"mis\") > 0.99'",
					"-O z", "-o", tmp_vcf)
command2    = paste("bcftools", "index", "--tbi", "-f", tmp_vcf)

system(command1)
system(command2)

tmp_gttable = sub("vcf.gz", "txt", tmp_vcf)
vcf_to_text = paste("bcftools", "query", "-H", "-s", paste(wgs_list, collapse = ","), tmp_vcf, "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID\\t%AF[\\t%GT]\\n"', "-o", tmp_gttable)

system(vcf_to_text)

gttable  = read_gttable(tmp_gttable, wgs_list)
gt_info  = gttable[, c("chrom", "pos", "ref", "alt", "rsid", "id", "af")]
gtmatrix = as.matrix(gttable[, wgs_list])
gtmatrix = suppressMessages(mapvalues(gtmatrix, from=c("0/0", "0/1", "1/0", "1/1", "./."), to=c(0, 0.5, 0.5, 1, NA)))

if (nrow(gt_info) > 0)
{
	class(gtmatrix) = "numeric" 
	gtmatrix        = gtmatrix[unlist(apply(gtmatrix, 1, FUN = function(x){sd(x, na.rm = TRUE)})) > 0,]
	gttable         = as.data.frame(gtmatrix)
	gt_info$rsid    = unlist(lapply(gt_info$rsid, function(x)
	{
		out = "."
		if(grepl(";", x) == TRUE){out = unlist(strsplit(x, ";"))[[2]]}

		return(out)
	}))

	write.table(gt_info , paste(outfolder, paste("gt_info", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	write.table(gtmatrix, paste(outfolder, paste("gt_data", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = TRUE , col.names = NA  )
}

system(paste("rm", "-r", tmp_folder))

# Write isoforms
isoform_info = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
isofolder    = paste(getwd(), "pipeline/1.3.genotype/use_isoform", sep = "/")
isoforms     = isoform_info[isoform_info$gene_id == gene_id,]

invisible(lapply(isoforms$transcript_id, function(transcript_id)
{
	suppressWarnings(file.symlink(from = paste(outfolder, paste("gt_info", gene_id, "txt", sep = "."), sep = "/"), to = paste(isofolder, paste("gt_info", transcript_id, "txt", sep = "."), sep = "/")))
	suppressWarnings(file.symlink(from = paste(outfolder, paste("gt_data", gene_id, "txt", sep = "."), sep = "/"), to = paste(isofolder, paste("gt_data", transcript_id, "txt", sep = "."), sep = "/")))
}))
