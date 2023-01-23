
setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")
source("script/functions.R" )
source("script/cibersort.R" )

library(Seurat)

annots_facs    = fread("input/tabula_muris/annotations_facs.csv"   , sep = ",", header = TRUE, data.table = FALSE)
meta_facs      = fread("input/tabula_muris/metadata_facs.csv"      , sep = ",", header = TRUE, data.table = FALSE)

load("input/tabula_muris/facs_Heart_seurat_tiss.Robj")
heart = tiss

heart = UpdateSeuratObject(heart)

heart = RunUMAP(heart, verbose = FALSE, dims = 1:20)

heart@meta.data[is.na(heart@meta.data$cell_ontology_class) == TRUE, "cell_ontology_class"] = "cardiac neuron"

saveRDS(object = heart, file = "pipeline/2.1.scrna_seq/tabula_muris.heart.rds")

# Save heart metadata

umapdata           = as.data.frame(heart@reductions$umap@cell.embeddings)
colnames(umapdata) = paste("umap", 1:2, sep = "")
tsnedata           = as.data.frame(heart@reductions$tsne@cell.embeddings)
colnames(tsnedata) = paste("tsne", 1:2, sep = "")
scmeta             = heart@meta.data
umapdata$barcode   = rownames(umapdata)
tsnedata$barcode   = rownames(tsnedata)
scmeta$barcode     = rownames(scmeta)
scmeta             = merge(scmeta, umapdata)
scmeta             = merge(scmeta, tsnedata)

fwrite(scmeta, "pipeline/2.1.scrna_seq/tabula_muris.heart.metadata.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# Deconvolution

geneinfo =              fread("input/phenotypes/gene_info.txt"                , sep = "\t", header = TRUE, data.table = FALSE)
metadata =              fread("pipeline/1.1.metadata/meta_rna.txt"            , sep = "\t", header = TRUE , data.table = FALSE)
tpm      = add_rownames(fread("pipeline/1.2.expression/tpm_gene.expressed.txt", sep = "\t", header = TRUE , data.table = FALSE))
metadata = metadata[metadata$run %in% colnames(tpm),]

rownames(geneinfo) = geneinfo$gene_id
gene2expr          = geneinfo[geneinfo$gene_id %in% rownames(tpm),]
genes_unique       = sort(table(gene2expr$gene_name), decreasing = TRUE)
genes_unique       = names(genes_unique[genes_unique == 1])
gene2expr          = gene2expr[gene2expr$gene_name %in% genes_unique,]
expr               = tpm[gene2expr$gene_id,]
rownames(expr    ) = gene2expr$gene_name

run_cibersort = function(scdata, name, tissues, metadata, expr)
{
    message(name)
    Idents(scdata)        = "cell_ontology_class"
	markers               = FindAllMarkers(scdata, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.1)
	markers               = markers[markers$p_val_adj <= 0.1,]
	genes2cibersort       = unique(markers$gene)
	sig.avg_exp           = AverageExpression(scdata, assays = "RNA", features = unique(genes2cibersort))
	sig.avg_exp           = sig.avg_exp$RNA
	colnames(sig.avg_exp) = gsub(" ", "_", colnames(sig.avg_exp))
	markers$gene          = toupper(markers$gene)
	rownames(sig.avg_exp) = toupper(rownames(sig.avg_exp))
	signature_gene_matrix = sig.avg_exp
	this_metadata         = metadata[metadata$tissue %in% tissues,]
	gene_ids              = intersect(rownames(signature_gene_matrix), rownames(expr))
	mixture_matrix        = expr                 [gene_ids,this_metadata$run]
	signature_matrix      = signature_gene_matrix[gene_ids,]
	cibersort             = CIBERSORT(signature_matrix, mixture_matrix, QN = FALSE, perm = 0)
	cibersort             = as.data.frame(cibersort)
	
    fwrite(signature_gene_matrix, paste("pipeline/2.1.scrna_seq/signature_gene_average_expression", name, "txt", sep = "."), sep = "\t", col.names = TRUE, row.names = TRUE)
    fwrite(markers              , paste("pipeline/2.1.scrna_seq/markers"                          , name, "txt", sep = "."), sep = "\t", col.names = TRUE, row.names = FALSE)
    fwrite(cibersort            , paste("pipeline/2.1.scrna_seq/cibersort"                        , name, "txt", sep = "."), sep = "\t", col.names = TRUE, row.names = TRUE )
}

#run_cibersort(heart, "heart"  , c("heart", "ipsc_cvpc"), metadata, expr)
#run_cibersort(heart, "arteria", c("arteria"           ), metadata, expr)

# Reanalyze cibersort and combine

cibersort_arteria = add_rownames(fread("pipeline/2.1.scrna_seq//cibersort.arteria.txt", sep = "\t", header = TRUE , data.table = FALSE))
cibersort_heart   = add_rownames(fread("pipeline/2.1.scrna_seq//cibersort.heart.txt"  , sep = "\t", header = TRUE , data.table = FALSE))

cibersort_arteria[,c("P-value", "Correlation", "RMSE")] = NULL
cibersort_heart  [,c("P-value", "Correlation", "RMSE")] = NULL

cell_types = sort(unique(c(colnames(cibersort_arteria), colnames(cibersort_heart))))

cibersort_heart  [setdiff(cell_types, colnames(cibersort_heart  ))] = 0
cibersort_arteria[setdiff(cell_types, colnames(cibersort_arteria))] = 0

cibersort = rbind(cibersort_heart[,cell_types], cibersort_arteria[,cell_types])
cibersort = cibersort[,order(colMeans(cibersort), decreasing = TRUE)]
colnames(cibersort)[colnames(cibersort) == "leukocyte"] = "immune_cell"

colnames(cibersort) = gsub(" ", "_", colnames(cibersort))
cibersort_combined  = cibersort
colnames(cibersort) = paste("regular", colnames(cibersort), sep = ".")

cibersort_combined$combined.cardiac_muscle_cell = cibersort_combined$cardiac_muscle_cell
cibersort_combined$combined.smooth_muscle_cell  = cibersort_combined$smooth_muscle_cell 
cibersort_combined$combined.immune_cell         = cibersort_combined$immune_cell 
cibersort_combined$combined.cardiac_neuron      = cibersort_combined$cardiac_neuron      
cibersort_combined$combined.endothelial_cell    = rowSums(cibersort_combined[,c("endocardial_cell", "endothelial_cell"  )])
cibersort_combined$combined.fibroblast          = rowSums(cibersort_combined[,c("fibroblast"      , "myofibroblast_cell")])

cibersort_combined = cibersort_combined[,grepl("combined", colnames(cibersort_combined)) == TRUE]
cibersort_combined = cibersort_combined[,order(colMeans(cibersort_combined), decreasing = TRUE)]

colnames(cibersort         ) = gsub("_cell", "", colnames(cibersort         ))
colnames(cibersort_combined) = gsub("_cell", "", colnames(cibersort_combined))

cibersort_all = cbind(cibersort, cibersort_combined)

fwrite(cibersort_all, "pipeline/2.1.scrna_seq/cibersort_results.txt", sep = "\t", col.names = TRUE, row.names = TRUE )

