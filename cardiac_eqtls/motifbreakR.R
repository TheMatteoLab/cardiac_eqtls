library(data.table)
library(ggplot2)
library(gplots)
source('~/scripts/functions.R')
library(grid)
library(ggrepel)
library(dplyr)

setwd("/projects/CARDIPS/analysis/cardiac_eqtls")
results.dir = "~/projects/cardiac_eqtls/TF_analyses/"

library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19) 

######### 1. Prepping Data ---------------------------------

# Input
eqtls.cvpc_adult = fread(paste0(results.dir, "homer/bed4fasta/eqtls.finemap.cvpc_adult.no_exons.txt"))

message(paste("# starting snps:", length(eqtls.cvpc_adult$rsid), length(unique(eqtls.cvpc_adult$rsid))))

# Remove variants with no rsids 
message("removing variants with no rsids....")
eqtls.cvpc_adult = eqtls.cvpc_adult[grepl("rs", eqtls.cvpc_adult$rsid),]

message(paste("# snps:", length(eqtls.cvpc_adult$rsid), length(unique(eqtls.cvpc_adult$rsid))))

# Remove indels
message("removing indels...")
eqtls.cvpc_adult = eqtls.cvpc_adult[nchar(eqtls.cvpc_adult$ref) == 1 & nchar(eqtls.cvpc_adult$alt) == 1,]

message(paste("# snps:", length(eqtls.cvpc_adult$rsid), length(unique(eqtls.cvpc_adult$rsid))))

message(paste("# input snps:", length(eqtls.cvpc_adult$rsid), length(unique(eqtls.cvpc_adult$rsid))))

######### 2. Running motifbreakR  ---------------------------------

# Get Genomic Ranges from rsids
snps.mb <- snps.from.rsid(rsid = eqtls.cvpc_adult$rsid,
                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)

message(paste("# snps.mb:", length(unique(names(snps.mb)))))

message("snps in snp.mb but not input")
names(snps.mb)[which(!names(snps.mb) %in% eqtls.cvpc_adult$rsid)]

# Running motifbreakR

message("running motifbreakR")
data(hocomoco)

results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

save(results, file = paste0(results.dir, "motifbreakR/eqtls.cvpc_adult.results"))

message("calculating pvalues...")
#load(paste0(results.dir, "motifbreakR/eqtls.cvpc_adult.finemap.motifbreakR.out"))

out = data.frame(results)

datalist = list()

for (rsid in unique(out$SNP_id)){
    	message(rsid)
	id = results[names(results) %in% rsid,]
    	out = tryCatch(
		{ data.frame(calculatePvalue(id))
		},
		error = function(cond) {
		message(paste("error:", rsid))
		message(cond)
		return(NA) 
		},
		warning = function(cond) {
		message(paste("warning:", rsid))
		message(cond)
		return(NULL)
		},
		finally = {
		message(paste(rsid))
		}	
	)	
	datalist[[rsid]] = out
}

pvalues = as.data.frame(rbindlist(datalist), stringsAsFactors = FALSE)

save(pvalues, file = paste0(results.dir, "motifbreakR/eqtls.cvpc_adult.finemap.motifbreakR.pvalues.out"))

print("finished!")


