
setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )


suppressPackageStartupMessages(library(stringr))

cms   = data.frame(udid = c("UDID011", "UDID033", "UDID097", "UDID122", "UDID126", "UDID134", "UDID135", "UDID154", "UDID157", "UDID171"), cell = "cm"  )
epdcs = data.frame(udid = c("UDID036", "UDID046", "UDID060", "UDID069", "UDID089", "UDID132", "UDID155", "UDID188", "UDID203", "UDID279"), cell = "epdc")
udids = rbind(cms, epdcs)

udids$atac_file = paste("input/epigenome/atac", udids$udid, "bam", sep = ".")

fwrite(udids, "pipeline/3.3.ase_test/udids.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

eqtl_genes    = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.gene.egenes.txt"         , sep = "\t", header = TRUE, data.table = FALSE)
int_genes     = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.gene.interactions.txt"   , sep = "\t", header = TRUE, data.table = FALSE)
fm_genes      = fread("pipeline/3.2.eqtls/eqtls_fine_map//cardiac_eqtls.gene.txt"      , sep = "\t", header = TRUE, data.table = FALSE)
eqtl_isoforms = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.isoform.egenes.txt"      , sep = "\t", header = TRUE, data.table = FALSE)
int_isoforms  = fread("pipeline/3.2.eqtls/eqtls/cardiac_eqtls.isoform.interactions.txt", sep = "\t", header = TRUE, data.table = FALSE)
fm_isoforms   = fread("pipeline/3.2.eqtls/eqtls_fine_map//cardiac_eqtls.isoform.txt"   , sep = "\t", header = TRUE, data.table = FALSE)

# select only cardiac muscle and the same number of variants with pval == 1
int_genes_test    = int_genes   [int_genes   $cell_specific %in% paste("cibersort.regular", c("cardiac_muscle"), sep = "."),]
int_isoforms_test = int_isoforms[int_isoforms$cell_specific %in% paste("cibersort.regular", c("cardiac_muscle"), sep = "."),]

int_genes_bg      = int_genes   [int_genes   $cell_specific == "" & int_genes   $pval_int > 0.5 & !int_genes   $transcript_id %in% int_genes_test   $transcript_id,]
int_isoforms_bg   = int_isoforms[int_isoforms$cell_specific == "" & int_isoforms$pval_int > 0.5 & !int_isoforms$transcript_id %in% int_isoforms_test$transcript_id,]

set.seed(1)
int_genes_bg      = int_genes_bg   [sample(1:nrow(int_genes_bg   ), size = nrow(int_genes_test   ), replace = FALSE),]
set.seed(2)
int_isoforms_bg   = int_isoforms_bg[sample(1:nrow(int_isoforms_bg), size = nrow(int_isoforms_test), replace = FALSE),]

int_genes_test   $test = TRUE
int_isoforms_test$test = TRUE
int_genes_bg     $test = FALSE
int_isoforms_bg  $test = FALSE

int_genes    = rbind(int_genes_test   , int_genes_bg   )
int_isoforms = rbind(int_isoforms_test, int_isoforms_bg)

eqtl_genes   $tr2type = paste(eqtl_genes   $transcript_id, eqtl_genes   $type)
int_genes    $tr2type = paste(int_genes    $transcript_id, int_genes    $type)
fm_genes     $tr2type = paste(fm_genes     $transcript_id, fm_genes     $type)
eqtl_isoforms$tr2type = paste(eqtl_isoforms$transcript_id, eqtl_isoforms$type)
int_isoforms $tr2type = paste(int_isoforms $transcript_id, int_isoforms $type)
fm_isoforms  $tr2type = paste(fm_isoforms  $transcript_id, fm_isoforms  $type)

# select only fine mapped variants from cell type specific associations
totest_genes    = merge(int_genes   [,c("tr2type", "transcript_id", "gene_id", "gene_name", "type", "chrom", "pos", "ref", "alt", "id", "af", "rsid", "beta", "se", "pval", "fdr", "cell_specific", "beta_int", "se_int", "pval_int", "qval_int", "test")], fm_genes   [,c("tr2type", "pos", "ref", "alt", "id", "af", "rsid", "pp")], by = "tr2type", suffixes = c("_lead", "_fm"))
totest_isoforms = merge(int_isoforms[,c("tr2type", "transcript_id", "gene_id", "gene_name", "type", "chrom", "pos", "ref", "alt", "id", "af", "rsid", "beta", "se", "pval", "fdr", "cell_specific", "beta_int", "se_int", "pval_int", "qval_int", "test")], fm_isoforms[,c("tr2type", "pos", "ref", "alt", "id", "af", "rsid", "pp")], by = "tr2type", suffixes = c("_lead", "_fm"))
totest_genes    = totest_genes   [nchar(totest_genes   $ref_fm) == 1 & nchar(totest_genes   $alt_fm) == 1,]
totest_isoforms = totest_isoforms[nchar(totest_isoforms$ref_fm) == 1 & nchar(totest_isoforms$alt_fm) == 1,]

rewrite_ase = function(ase, totest, udids)
{
    totest$coord = paste(paste("chr", totest$chrom, sep = ""), totest$pos_fm)
    ase$coord    = paste(ase$chrom   , ase$pos)
    ase          = merge(unique(totest[,c("transcript_id", "gene_id", "gene_name", "coord", "id_lead", "rsid_lead", "af_lead", "id_fm", "rsid_fm", "af_fm", "ref_fm", "alt_fm", "test")]), ase[,c("coord", "udid", "cov", "gts")])
	ase$gts      = gsub("\\*", "", ase$gts)
    
    if(nrow(ase) > 0)
    {
        ase$gts      = toupper(ase$gts)
        ase          = cbind(ase, as.data.frame(rbindlist(lapply(1:nrow(ase), function(ii)
        {
            myref = ase[ii, "ref_fm"]
            myalt = ase[ii, "alt_fm"]
            gts   = ase[ii, "gts"   ]
			
			#message(paste(gts, myref, myalt))
			
			if(nchar(gts) > 0)
			{
				if( myref %in% c("A", "C", "G", "T")){ref_n = str_count(gts, pattern = myref)}
				if( myalt %in% c("A", "C", "G", "T")){alt_n = str_count(gts, pattern = myalt)}
				if(!myref %in% c("A", "C", "G", "T")){ref_n = 0}
				if(!myalt %in% c("A", "C", "G", "T")){alt_n = 0}
			}else
			{
				ref_n = 0
				alt_n = 0
			}
            out   = data.frame(ref_n = ref_n, alt_n = alt_n)

            return(out)
        })), stringsAsFactors = FALSE))

        #ase = ase[ase$ref_n > 0 & ase$alt_n > 0,]
		
		#ase$major = unlist(lapply(1:nrow(ase), function(ii){max(as.numeric(ase[ii, c("ref_n", "alt_n")]))}))
		#ase$minor = unlist(lapply(1:nrow(ase), function(ii){min(as.numeric(ase[ii, c("ref_n", "alt_n")]))}))
		#ase$maf   = ase$minor / rowSums(ase[,c("major", "minor")])
		ase$af   = ase$alt_n / rowSums(ase[,c("ref_n", "alt_n")])

    }else
    {
        #ase = cbind(ase, data.frame(ref_n = numeric(), alt_n = numeric(), major = numeric(), minor = numeric(), maf = numeric()))
        ase = cbind(ase, data.frame(ref_n = numeric(), alt_n = numeric(), af = numeric()))
    }
	
	ase = merge(ase, udids[,c("udid", "cell")])
    
    return(ase)
}

mpileup = function(totest, name, udid, udids, assay = "atac")
{
    message(paste(name, udid))
    
    positions       = unique(totest[,c("chrom", "pos_fm")])
    positions$chrom = paste("chr", positions$chrom, sep = "")
    positions       = positions[order(positions$chrom, positions$pos_fm),]
    positions_file  = paste("pipeline/3.3.ase_test/tmp", paste("TMP", assay, name, udid, "txt"    , sep = "."), sep = "/")
    mpileup_file    = paste("pipeline/3.3.ase_test/tmp", paste("TMP", assay, name, udid, "mpileup", sep = "."), sep = "/")
    bam_file        = udids[udids$udid == udid, paste(assay, "file", sep = "_")]
    command         = paste("samtools", "mpileup", "--excl-flags", 3852, "-l", positions_file, "-a", "-o", mpileup_file, bam_file)
    
    fwrite(positions, positions_file, sep = "\t", col.names = FALSE, row.names = FALSE)
    
    #system(command)
    
    if(file.size(mpileup_file) > 0)
    {
        indata           = fread(mpileup_file, sep = "\t", header = FALSE, data.table = FALSE)[,c(1,2,4,5)]
        colnames(indata) = c("chrom", "pos", "cov", "gts")
        indata$udid      = udid
    }else
    {
        indata = data.frame(chrom = numeric(), pos = numeric(), cov = numeric(), gts = character(), udid = character()) 
    }
    
    fwrite(indata , paste("pipeline/3.3.ase_test/udids", paste(assay, name, udid, "ase_by_udid", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
    return(indata)
}

ase_genes    = rewrite_ase(as.data.frame(rbindlist(lapply(udids$udid, function(udid){mpileup(totest_genes   , "gene"   , udid, udids, "atac")})), stringsAsFactors = FALSE), totest_genes   , udids)
#ase_isoforms = rewrite_ase(as.data.frame(rbindlist(lapply(udids$udid, function(udid){mpileup(totest_isoforms, "isoform", udid, udids, "atac")})), stringsAsFactors = FALSE), totest_isoforms, udids)

fwrite(ase_genes   , paste("pipeline/3.3.ase_test", paste("atac", "gene"   ,"ase_by_udid", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
#fwrite(ase_isoforms, paste("pipeline/3.3.ase_test", paste("atac", "isoform","ase_by_udid", "txt", sep = "."), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE)
