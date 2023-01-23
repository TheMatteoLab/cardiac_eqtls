
setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )


gene_info               = fread("input/phenotypes/gene_info.txt"   , sep = "\t", header = TRUE, data.table = FALSE)
isof_info               = fread("input/phenotypes/isoform_info.txt", sep = "\t", header = TRUE, data.table = FALSE)
gene_info$transcript_id = gene_info$gene_id

diffgene = fread(paste("pipeline/4.1.differential_expression", "diffexp.txt"     , sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
diffcell = fread(paste("pipeline/4.1.differential_expression", "diffexp_cell.txt", sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)

diffgene$diffexp = FALSE
diffgene[diffgene$qval < 0.05, "diffexp"] = TRUE

gmt2list = function(gs)
{    
    indata = readLines(paste("input/phenotypes/msigdb", paste(gs, "v7.1.symbols.gmt", sep = "."), sep = "/"))
    indata = lapply(indata, function(x)
    {
        this = unlist(strsplit(x, "\t"))
        out  = list(gene_set = this[[1]], url = this[[2]], gene_names = this[3:length(this)])
        
        return(out)
    })
    names(indata) = unlist(lapply(indata, function(x){x$gene_set}))
    
    return(indata)
}

msigdb               = c("c2.cp.biocarta", "c2.cp.kegg", "c2.cp.reactome", "c5.bp", "c5.cc", "c5.mf", "h.all")
#msigdb               = c("h.all")
genesets2test        = lapply(msigdb, gmt2list)
names(genesets2test) = msigdb

calculate_functional_enrichment_tissue_by_gene_set = function(tissue1, tissue2, type, gs, gene_set, diffgene)
{
    diffexp    = diffgene[diffgene$type == type & diffgene$tissue1 == tissue1 & diffgene$tissue2 == tissue2,]
    out        = data.frame(tissue1 = tissue1, tissue2 = tissue2, type = type, gs_source = gs, gene_set = gene_set[["gene_set"]], url = gene_set[["url"]], ngenes = length(gene_set[["gene_names"]]))
    genes_in   = gene_set[["gene_names"]]
    diffexp$gs = FALSE
    
    if(nrow(diffexp[diffexp$gene_name %in% genes_in,]) > 1)
    {
        diffexp[diffexp$gene_name %in% genes_in, "gs"] = TRUE

        test = t.test(diffexp[diffexp$gs == TRUE, "beta"], diffexp[diffexp$gs == FALSE, "beta"])
        out  = cbind(out, data.frame(ttest_estimate_in = test$estimate[[1]], ttest_estimate_out = test$estimate[[2]], ttest_ci1 = test$conf.int[[1]], ttest_ci2 = test$conf.int[[2]], ttest_pval = test$p.value))

        tofisher = matrix(0, nrow = 3, ncol = 2)
        rownames(tofisher) = c("over", "under", "no")
        colnames(tofisher) = c("yes", "no")

        tofisher["over" , "yes"] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta > 0 & diffexp$gs == TRUE ,])
        tofisher["over" , "no" ] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta > 0 & diffexp$gs == FALSE,])
        tofisher["under", "yes"] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta < 0 & diffexp$gs == TRUE ,])
        tofisher["under", "no" ] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta < 0 & diffexp$gs == FALSE,])
        tofisher["no"   , "yes"] = nrow(diffexp[diffexp$qval > 0.05                    & diffexp$gs == TRUE ,])
        tofisher["no"   , "no" ] = nrow(diffexp[diffexp$qval > 0.05                    & diffexp$gs == FALSE,])

        test_up = fisher.test(tofisher[c("over" , "no"),])
        test_dn = fisher.test(tofisher[c("under", "no"),])

        out  = cbind(out, data.frame(fisher_up_estimate = test_up$estimate, fisher_up_ci1 = test_up$conf.int[[1]], fisher_up_ci2 = test_up$conf.int[[2]], fisher_up_pval = test_up$p.value))
        out  = cbind(out, data.frame(fisher_dn_estimate = test_dn$estimate, fisher_dn_ci1 = test_dn$conf.int[[1]], fisher_dn_ci2 = test_dn$conf.int[[2]], fisher_dn_pval = test_dn$p.value))
    }else
    {
        out  = cbind(out, data.frame(ttest_estimate_in  = NA, ttest_estimate_out = NA, ttest_ci1     = NA, ttest_ci2      = NA, ttest_pval = NA))
        out  = cbind(out, data.frame(fisher_up_estimate = NA, fisher_up_ci1      = NA, fisher_up_ci2 = NA, fisher_up_pval = NA))
        out  = cbind(out, data.frame(fisher_dn_estimate = NA, fisher_dn_ci1      = NA, fisher_dn_ci2 = NA, fisher_dn_pval = NA))

    }
    return(out)
}

#tests = as.data.frame(rbindlist(lapply(sort(unique(diffgene$type)), function(type)
tests = as.data.frame(rbindlist(lapply(c("gene_tpm", "isoform_use"), function(type)
{
    as.data.frame(rbindlist(lapply(msigdb, function(gs)
    {
        message(paste(type, gs))
        gene_sets = genesets2test[[gs]]
        as.data.frame(rbindlist(lapply(gene_sets, function(gene_set)
        {
            #message(gene_set[["gene_set"]])
			#calculate_functional_enrichment_tissue_by_gene_set("ipsc_cvpc", "heart"  , type, gs, gene_set, diffgene)
            out = rbind(calculate_functional_enrichment_tissue_by_gene_set("ipsc_cvpc", "heart"  , type, gs, gene_set, diffgene),
                        calculate_functional_enrichment_tissue_by_gene_set("ipsc_cvpc", "arteria", type, gs, gene_set, diffgene),
                        calculate_functional_enrichment_tissue_by_gene_set("heart"    , "arteria", type, gs, gene_set, diffgene)
                       )
            return(out)
        })), stringsAsFactors = FALSE)
    })), stringsAsFactors = FALSE)
})), stringsAsFactors = FALSE)

tests = tests[is.na(tests$ttest_estimate_in) == FALSE, ]

tests$ttest_fdr     = p.adjust(tests$ttest_pval    , method = "BH")
tests$fisher_up_fdr = p.adjust(tests$fisher_up_pval, method = "BH")
tests$fisher_dn_fdr = p.adjust(tests$fisher_dn_pval, method = "BH")

fwrite(tests, "pipeline/4.1.differential_expression/functional_enrichment.tissue.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

tests_tissue = tests

calculate_functional_enrichment_cell_by_gene_set = function(cell_type, type, gs, gene_set, diffcell)
{
    diffexp    = diffcell[diffcell$type == type & diffcell$cell_type == cell_type,]
    out        = data.frame(cell_type = cell_type, type = type, gs_source = gs, gene_set = gene_set[["gene_set"]], url = gene_set[["url"]], ngenes = length(gene_set[["gene_names"]]))
    genes_in   = gene_set[["gene_names"]]
    diffexp$gs = FALSE
    
    if(nrow(diffexp[diffexp$gene_name %in% genes_in,]) > 1)
    {
        diffexp[diffexp$gene_name %in% genes_in, "gs"] = TRUE

        test = t.test(diffexp[diffexp$gs == TRUE, "beta"], diffexp[diffexp$gs == FALSE, "beta"])
        out  = cbind(out, data.frame(ttest_estimate_in = test$estimate[[1]], ttest_estimate_out = test$estimate[[2]], ttest_ci1 = test$conf.int[[1]], ttest_ci2 = test$conf.int[[2]], ttest_pval = test$p.value))

        tofisher = matrix(0, nrow = 3, ncol = 2)
        rownames(tofisher) = c("over", "under", "no")
        colnames(tofisher) = c("yes", "no")

        tofisher["over" , "yes"] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta > 0 & diffexp$gs == TRUE ,])
        tofisher["over" , "no" ] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta > 0 & diffexp$gs == FALSE,])
        tofisher["under", "yes"] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta < 0 & diffexp$gs == TRUE ,])
        tofisher["under", "no" ] = nrow(diffexp[diffexp$qval < 0.05 & diffexp$beta < 0 & diffexp$gs == FALSE,])
        tofisher["no"   , "yes"] = nrow(diffexp[diffexp$qval > 0.05                    & diffexp$gs == TRUE ,])
        tofisher["no"   , "no" ] = nrow(diffexp[diffexp$qval > 0.05                    & diffexp$gs == FALSE,])

        test_up = fisher.test(tofisher[c("over" , "no"),])
        test_dn = fisher.test(tofisher[c("under", "no"),])

        out  = cbind(out, data.frame(fisher_up_estimate = test_up$estimate, fisher_up_ci1 = test_up$conf.int[[1]], fisher_up_ci2 = test_up$conf.int[[2]], fisher_up_pval = test_up$p.value))
        out  = cbind(out, data.frame(fisher_dn_estimate = test_dn$estimate, fisher_dn_ci1 = test_dn$conf.int[[1]], fisher_dn_ci2 = test_dn$conf.int[[2]], fisher_dn_pval = test_dn$p.value))
    }else
    {
        out  = cbind(out, data.frame(ttest_estimate_in  = NA, ttest_estimate_out = NA, ttest_ci1     = NA, ttest_ci2      = NA, ttest_pval = NA))
        out  = cbind(out, data.frame(fisher_up_estimate = NA, fisher_up_ci1      = NA, fisher_up_ci2 = NA, fisher_up_pval = NA))
        out  = cbind(out, data.frame(fisher_dn_estimate = NA, fisher_dn_ci1      = NA, fisher_dn_ci2 = NA, fisher_dn_pval = NA))
    }
    return(out)
}

cell_types = sort(unique(diffcell$cell_type))
cell_types = cell_types[grepl("regular", cell_types) == TRUE]

#tests = as.data.frame(rbindlist(lapply(sort(unique(diffgene$type)), function(type)
tests = as.data.frame(rbindlist(lapply(c("gene_tpm", "isoform_use"), function(type)
{
    as.data.frame(rbindlist(lapply(cell_types, function(cell_type)
    {
        as.data.frame(rbindlist(lapply(msigdb, function(gs)
        {
            message(paste(type, cell_type, gs))
            gene_sets = genesets2test[[gs]]
            as.data.frame(rbindlist(lapply(gene_sets, function(gene_set)
            {
                #message(gene_set[["gene_set"]])
                calculate_functional_enrichment_cell_by_gene_set(cell_type, type, gs, gene_set, diffcell)
            })), stringsAsFactors = FALSE)
        })), stringsAsFactors = FALSE)
    })), stringsAsFactors = FALSE)
})), stringsAsFactors = FALSE)

tests = tests[is.na(tests$ttest_estimate_in) == FALSE, ]

tests$ttest_fdr     = p.adjust(tests$ttest_pval    , method = "BH")
tests$fisher_up_fdr = p.adjust(tests$fisher_up_pval, method = "BH")
tests$fisher_dn_fdr = p.adjust(tests$fisher_dn_pval, method = "BH")

fwrite(tests, "pipeline/4.1.differential_expression/functional_enrichment.cell.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
