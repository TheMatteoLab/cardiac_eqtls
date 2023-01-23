setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )

dir.create("pipeline/1.2.expression"            , showWarnings = FALSE)
dir.create("pipeline/1.2.expression/tpm_gene"   , showWarnings = FALSE)
dir.create("pipeline/1.2.expression/tpm_isoform", showWarnings = FALSE)
dir.create("pipeline/1.2.expression/use_isoform", showWarnings = FALSE)

# input files

metadata_original                    =              fread("pipeline/1.1.metadata/metadata.txt"              , sep = "\t", header = TRUE, data.table = FALSE)
meta_rna_original                    =              fread("pipeline/1.1.metadata/meta_rna.txt"              , sep = "\t", header = TRUE, data.table = FALSE)
meta_subject_original                =              fread("pipeline/1.1.metadata/meta_subject.txt"          , sep = "\t", header = TRUE, data.table = FALSE)
flagstat_ipscore                     = add_rownames(fread("input/phenotypes/ipscore.flagstat.txt"           , sep = "\t", header = TRUE, data.table = FALSE))
gene_tpm_ipscore                     = add_rownames(fread("input/phenotypes/ipscore.gene_tpm.txt"           , sep = "\t", header = TRUE, data.table = FALSE))
isoform_percent_use_ipscore          = add_rownames(fread("input/phenotypes/ipscore.isoform_percent_use.txt", sep = "\t", header = TRUE, data.table = FALSE))
isoform_tpm_ipscore                  = add_rownames(fread("input/phenotypes/ipscore.isoform_tpm.txt"        , sep = "\t", header = TRUE, data.table = FALSE))
flagstat_gtex                        = add_rownames(fread("input/phenotypes/gtex.flagstat.txt"              , sep = "\t", header = TRUE, data.table = FALSE))
gene_tpm_gtex                        = add_rownames(fread("input/phenotypes/gtex.gene_tpm.txt"              , sep = "\t", header = TRUE, data.table = FALSE))
isoform_percent_use_gtex             = add_rownames(fread("input/phenotypes/gtex.isoform_percent_use.txt"   , sep = "\t", header = TRUE, data.table = FALSE))
isoform_tpm_gtex                     = add_rownames(fread("input/phenotypes/gtex.isoform_tpm.txt"           , sep = "\t", header = TRUE, data.table = FALSE))

## Adjust this part when we have all GTEx

flagstat            = rbind(flagstat_ipscore           , flagstat_gtex           )
gene_tpm            = cbind(gene_tpm_ipscore           , gene_tpm_gtex           )
isoform_percent_use = cbind(isoform_percent_use_ipscore, isoform_percent_use_gtex)
isoform_tpm         = cbind(isoform_tpm_ipscore        , isoform_tpm_gtex        )

ids                 = colnames(gene_tpm)

# Get expression matrices and filter metadata files

meta_rna     = meta_rna_original    [meta_rna_original$run            %in% ids                        ,]
metadata     = metadata_original    [metadata_original$run            %in% ids                        ,]
meta_subject = meta_subject_original[meta_subject_original$subject_id %in% unique(metadata$subject_id),]

flagstat                = flagstat           [ meta_rna$run,]
gene_tpm                = gene_tpm           [,meta_rna$run ]
isoform_percent_use     = isoform_percent_use[,meta_rna$run ]
isoform_tpm             = isoform_tpm        [,meta_rna$run ]

# Create RNA covariates (from flagstat)
covariates_rna = data.frame(row.names   = rownames(flagstat),
							total_reads = flagstat$total_reads,
							uniquely_mapped_reads                         = 100 * flagstat$both_mapped / flagstat$total_reads,
							uniquely_mapped_reads_to_canonical_chromsomes = 100 * rowSums(flagstat[, c("autosomal_reads", "chrX_reads", "chrY_reads")]) / flagstat$total_reads,
							mitochondrial_reads                           = 100 * flagstat$chrM_reads / flagstat$total_reads)


# Find expressed genes
gene_info    = fread("input/phenotypes/gene_info.txt"   , sep = "\t", header = TRUE, data.table = FALSE)
isoform_info = fread("input/phenotypes/isoform_info.txt", sep = "\t", header = TRUE, data.table = FALSE)

gene_info$transcript_id = gene_info$gene_id

# Blacklist genes and isoforms: remove all genes with either gene expression or isoform use significantly associated with read length
bl_gene    = fread("/projects/CARDIPS/analysis/trimread_validation/pval/gene_tpm.txt"           , sep = " ", header = FALSE, data.table = FALSE)
bl_isoform = fread("/projects/CARDIPS/analysis/trimread_validation/pval/isoform_percent_use.txt", sep = " ", header = FALSE, data.table = FALSE)

expressed                                       = as.matrix(isoform_percent_use)
expressed[as.matrix(isoform_percent_use) <  10] = 0
expressed[as.matrix(isoform_percent_use) >= 10] = 1
tofilter                                        = data.frame(transcript_id = rownames(isoform_percent_use), use = rowSums(expressed), pass = FALSE)
tofilter[tofilter$use >= (0.1 * ncol(isoform_percent_use)), "pass"] = TRUE

blacklist_gene    = unique(  bl_gene[bl_gene$V3 <= 0.05, "V1"])
blacklist_isoform = unique(c(bl_gene[bl_gene$V3 <= 0.05, "V1"], isoform_info[isoform_info$transcript_id %in% bl_isoform[bl_isoform$V4 <= 0.05 & !bl_isoform$V2 %in% tofilter[tofilter$pass == FALSE, "transcript_id"], "V2"], "gene_id"]))

writeLines(blacklist_gene   , con  = "pipeline/1.2.expression/blacklist_gene.txt"   , sep = "\n")
writeLines(blacklist_isoform, con  = "pipeline/1.2.expression/blacklist_isoform.txt", sep = "\n")

gene_info    = gene_info   [!gene_info   $gene_id %in% blacklist_gene   , ]
isoform_info = isoform_info[!isoform_info$gene_id %in% blacklist_isoform, ]

gene_tpm            = gene_tpm           [gene_info   $gene_id      ,]
isoform_tpm         = isoform_tpm        [isoform_info$transcript_id,]
isoform_percent_use = isoform_percent_use[isoform_info$transcript_id,]

divide_phenotypes_tissue("pipeline/1.2.expression/tpm_gene", gene_tpm, geneinfo = gene_info, run_peer = TRUE, normalize = TRUE, gene_ids = NULL, filter_exp = TRUE, phenotype_min_value = 1, phenotype_min_samples = 0.1, peer_factor_n = 300, n_genes = 10000, use = FALSE)

#divide_phenotypes_tissue = function(outfolder, tpm, geneinfo, run_peer = FALSE, normalize = FALSE, gene_ids = NULL, filter_exp = FALSE, phenotype_min_value = 1, phenotype_min_samples = 0.1, peer_factor_n = 10, n_genes = 1000, use = FALSE)
#{
#    message("---------------------------")
#    
#    expdata = normalize_tpm(outfolder, tpm, geneinfo, normalize, gene_ids, filter_exp, min_tpm = phenotype_min_value, min_samples = phenotype_min_samples, use)
#    
#    if(run_peer == TRUE){peerdata = calculate_peer_factors(outfolder, expdata, peer_factor_n, n_genes)}
#    
#    #divide_phenotypes_by_gene(expdata, outfolder)
#}
#divide_phenotypes_tissue("pipeline/1.2.expression/tpm_gene", gene_tpm, geneinfo = gene_info, run_peer = FALSE, normalize = TRUE, gene_ids = NULL, filter_exp = TRUE, phenotype_min_value = 1, phenotype_min_samples = 0.1, peer_factor_n = 300, n_genes = 10000, use = FALSE)

genes_expressed        = readLines("pipeline/1.2.expression/tpm_gene.gene_ids.txt")
isoforms_expressed     = isoform_info[isoform_info$gene_id %in% genes_expressed,]
gene2isoform_expressed = table(isoforms_expressed$gene_id)
isoforms_to_remove     = names(gene2isoform_expressed[gene2isoform_expressed == 1])
isoforms_expressed     = isoforms_expressed[!isoforms_expressed$gene_id %in% isoforms_to_remove, "transcript_id"]


divide_phenotypes_tissue("pipeline/1.2.expression/tpm_isoform", isoform_tpm, geneinfo = isoform_info, run_peer = FALSE, normalize = TRUE, gene_ids = isoforms_expressed, filter_exp = TRUE, phenotype_min_value = 1, phenotype_min_samples = 0.1, use = FALSE)

isoforms_expressed = readLines("pipeline/1.2.expression/tpm_isoform.gene_ids.txt")

divide_phenotypes_tissue("pipeline/1.2.expression/use_isoform", isoform_percent_use, geneinfo = isoform_info, run_peer = FALSE, normalize = TRUE, gene_ids = isoforms_expressed, filter_exp = TRUE, phenotype_min_value = 10, phenotype_min_samples = 0.1, use = TRUE)

isoforms_use_expressed = readLines("pipeline/1.2.expression/use_isoform.gene_ids.txt")

# Rewrite gene and isoform information
gene_info_filtered    = gene_info   [gene_info   $gene_id       %in% genes_expressed       ,]
isoform_info_filtered = isoform_info[isoform_info$transcript_id %in% isoforms_use_expressed,]

fwrite(gene_info_filtered   , "pipeline/1.2.expression/gene_info.txt"    , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(isoform_info_filtered, "pipeline/1.2.expression/isoform_info.txt" , sep = "\t", col.names = TRUE, row.names = FALSE)

# Write metadata output
peer_factors   = add_rownames(fread("pipeline/1.2.expression/tpm_gene.peer.factors.txt", sep = "\t", header = TRUE , data.table = FALSE))
covariates_rna = merge(peer_factors, covariates_rna, by.x = "assay_id", by.y = "row.names")

wgs_id_vcf_ipscore = unlist(strsplit(tail(system("bcftools view -h input/genotypes/ipscore.vcf.gz", intern = TRUE), n = 1), split = "\t"))
wgs_id_vcf_gtex    = unlist(strsplit(tail(system("bcftools view -h input/genotypes/gtex.vcf.gz"   , intern = TRUE), n = 1), split = "\t"))
wgs_id_vcf         = c(wgs_id_vcf_ipscore[10:length(wgs_id_vcf_ipscore)], wgs_id_vcf_gtex[10:length(wgs_id_vcf_gtex)])
metadata$vcf       = FALSE

metadata[metadata$wgs_id %in% wgs_id_vcf, "vcf"] = TRUE

fwrite(covariates_rna, "pipeline/1.2.expression/covariates_rna.txt" , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(metadata      , "pipeline/1.2.expression/metadata.txt"       , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(meta_subject  , "pipeline/1.2.expression/meta_subject.txt"   , sep = "\t", col.names = TRUE, row.names = FALSE)

