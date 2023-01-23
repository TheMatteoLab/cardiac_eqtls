setwd("/frazer01/projects/CARDIPS/analysis/cardiac_eqtls")

source("script/packages.R"  )
source("script/input_data.R")

dir.create("pipeline/1.1.metadata", showWarnings = FALSE)


# input files
meta_subject_ipscore_original = fread("input/metadata/ipscore.subject_metadata.txt", sep = "\t", header = TRUE, data.table = FALSE)
meta_wgs_ipscore_original     = fread("input/metadata/ipscore.wgs_metadata.txt"    , sep = "\t", header = TRUE, data.table = FALSE)
meta_rna1_ipscore_original    = fread("input/metadata/ipscore.sample_metadata1.txt", sep = "\t", header = TRUE, data.table = FALSE)
meta_rna2_ipscore_original    = fread("input/metadata/ipscore.sample_metadata2.txt", sep = "\t", header = TRUE, data.table = FALSE)
meta_subject_gtex_original    = fread("input/metadata/gtex.subject_metadata.txt"   , sep = "\t", header = TRUE, data.table = FALSE, skip = 10)
meta_rna_gtex_original        = fread("input/metadata/gtex.sample_metadata.txt"    , sep = "\t", header = TRUE, data.table = FALSE)

addDashUUID = function(uuid)
{
    splitted = unlist(strsplit(uuid, ""))

    out = paste(paste(splitted[ 1: 8], collapse = ""),
                paste(splitted[ 9:12], collapse = ""),
                paste(splitted[13:16], collapse = ""),
                paste(splitted[17:20], collapse = ""),
                paste(splitted[21:32], collapse = ""),
                sep = "-"
               )
    return(out)
}

meta_wgs_ipscore_original$subject_id = unlist(lapply(meta_wgs_ipscore_original$subject_id, addDashUUID))
meta_wgs_ipscore_original$wgs_id     = unlist(lapply(meta_wgs_ipscore_original$id        , addDashUUID))

# combine iPSCORE and GTEx: RNA metadata
meta_diff_ipscore           = meta_rna1_ipscore_original[,c("UDID", "ipscore_id", "subject_uuid", "%cTnT")]
colnames(meta_diff_ipscore) = c("udid", "subject_name", "subject_id", "ctnt")
meta_rna_ipscore            = merge(meta_diff_ipscore, meta_rna2_ipscore_original[is.na(meta_rna2_ipscore_original$UDID) == FALSE, c("rna_assay_uuid", "UDID")], by.x = "udid", by.y = "UDID")
meta_rna_ipscore$assay_id   = meta_rna_ipscore$rna_assay_uuid
meta_rna_ipscore$run        = meta_rna_ipscore$rna_assay_uuid
meta_rna_ipscore$body_site  = "ipsc_cvpc"
meta_rna_ipscore$tissue     = "ipsc_cvpc"
meta_rna_ipscore$is_paired  = TRUE
meta_rna_ipscore$study      = "ipscore"
meta_rna_ipscore            = meta_rna_ipscore[,c("subject_id", "subject_name", "assay_id", "run", "tissue", "body_site", "study", "is_paired", "udid", "ctnt")]
meta_rna_gtex               = data.frame(subject_id       = meta_rna_gtex_original$submitted_subject_id,
										 subject_name     = meta_rna_gtex_original$submitted_subject_id,
										 assay_id         = meta_rna_gtex_original$biospecimen_repository_sample_id, 
										 run              = meta_rna_gtex_original$Run, 
										 tissue           = meta_rna_gtex_original$body_site, 
										 body_site        = meta_rna_gtex_original$body_site, 
										 study            = "gtex", 
										 is_paired        = TRUE, 
										 udid             = NA,
										 ctnt             = NA,
										 stringsAsFactors = FALSE)

meta_rna_gtex = meta_rna_gtex[meta_rna_gtex$tissue != "Pancreas",]
meta_rna_gtex[meta_rna_gtex$run %in% meta_rna_gtex_original[meta_rna_gtex_original$LibraryLayout == "SINGLE", "Run"], "is_paired"] = FALSE
meta_rna_gtex[meta_rna_gtex$body_site == "Heart - Left Ventricle"  , "body_site"] = "heart_ventricle"
meta_rna_gtex[meta_rna_gtex$body_site == "Heart - Atrial Appendage", "body_site"] = "heart_atrium"
meta_rna_gtex[meta_rna_gtex$body_site == "Artery - Aorta"          , "body_site"] = "arteria_aorta"
meta_rna_gtex[meta_rna_gtex$body_site == "Artery - Coronary"       , "body_site"] = "arteria_coronary"
meta_rna_gtex[meta_rna_gtex$tissue    == "Heart - Left Ventricle"  , "tissue"   ] = "heart"
meta_rna_gtex[meta_rna_gtex$tissue    == "Heart - Atrial Appendage", "tissue"   ] = "heart"
meta_rna_gtex[meta_rna_gtex$tissue    == "Artery - Aorta"          , "tissue"   ] = "arteria"
meta_rna_gtex[meta_rna_gtex$tissue    == "Artery - Coronary"       , "tissue"   ] = "arteria"

meta_rna_unfiltered = rbind(meta_rna_ipscore, meta_rna_gtex)

# combine iPSCORE and GTEx: subject and WGS metadata
meta_subject_gtex              = meta_subject_gtex_original[,c("SUBJID", "SEX", "AGE", "HGHT", "WGHT", "BMI")]
colnames(meta_subject_gtex)    = c("subject_id", "sex", "age", "height", "weight", "bmi")
meta_subject_gtex$wgs_id       = meta_subject_gtex$subject_id
meta_subject_gtex$subject_name = meta_subject_gtex$subject_id
meta_subject_gtex$sex          = suppressMessages(mapvalues(meta_subject_gtex$sex , from=c(1,2), to = c("M", "F")))

meta_subject_ipscore              = meta_subject_ipscore_original
meta_subject_ipscore$subject_name = gsub("CARDiPS", "iPSCORE", meta_subject_ipscore$ext_name)
meta_subject_ipscore              = merge(unique(meta_rna_ipscore[,c("subject_id", "subject_name")]), meta_subject_ipscore[,c("subject_name", "sex", "age", "height", "weight")])
meta_subject_ipscore$age          = as.integer(meta_subject_ipscore$age   )
meta_subject_ipscore$height       = as.numeric(meta_subject_ipscore$height)
meta_subject_ipscore$weight       = as.numeric(meta_subject_ipscore$weight)

meta_subject_ipscore[meta_subject_ipscore$sex == "M" & is.na(meta_subject_ipscore$height) == TRUE, "height"] = mean(meta_subject_ipscore[meta_subject_ipscore$sex == "M", "height"], na.rm = TRUE)
meta_subject_ipscore[meta_subject_ipscore$sex == "F" & is.na(meta_subject_ipscore$height) == TRUE, "height"] = mean(meta_subject_ipscore[meta_subject_ipscore$sex == "F", "height"], na.rm = TRUE)
meta_subject_ipscore[meta_subject_ipscore$sex == "M" & is.na(meta_subject_ipscore$weight) == TRUE, "weight"] = mean(meta_subject_ipscore[meta_subject_ipscore$sex == "M", "weight"], na.rm = TRUE)
meta_subject_ipscore[meta_subject_ipscore$sex == "F" & is.na(meta_subject_ipscore$weight) == TRUE, "weight"] = mean(meta_subject_ipscore[meta_subject_ipscore$sex == "F", "weight"], na.rm = TRUE)

meta_subject_ipscore$bmi = meta_subject_ipscore$weight / (meta_subject_ipscore$height)^2 * 703
meta_subject_ipscore     = merge(meta_subject_ipscore, meta_wgs_ipscore_original[meta_wgs_ipscore_original$cell != "iPSC" & meta_wgs_ipscore_original$status == 0, c("subject_id", "wgs_id")])

meta_subject_unfiltered = rbind(meta_subject_ipscore[,c("subject_id", "subject_name", "sex", "age", "height", "weight", "bmi")],
				                meta_subject_gtex   [,c("subject_id", "subject_name", "sex", "age", "height", "weight", "bmi")])


meta_wgs_unfiltered = unique(rbind(meta_subject_ipscore[,c("subject_id", "wgs_id")],
				                   meta_subject_gtex   [,c("subject_id", "wgs_id")]))

# Retain only samples that have RNA-seq data
assay_ids   = unique(meta_rna_unfiltered$assay_id  )
subject_ids = unique(intersect(meta_rna_unfiltered$subject_id, meta_subject_unfiltered$subject_id))

meta_rna         = meta_rna_unfiltered    [meta_rna_unfiltered    $subject_id %in% subject_ids,]
meta_subject     = meta_subject_unfiltered[meta_subject_unfiltered$subject_id %in% subject_ids,]
meta_wgs         = meta_wgs_unfiltered    [meta_wgs_unfiltered    $subject_id %in% subject_ids,]
metadata         = merge(meta_rna[,c("assay_id", "run", "subject_id", "study")], meta_wgs)
metadata         = metadata[,c("subject_id", "wgs_id", "assay_id", "run", "study")]
meta_subject$sex = suppressMessages(mapvalues(meta_subject$sex , from=c("M", "F"), to = c(0,1)))


fwrite(metadata    , "pipeline/1.1.metadata/metadata.txt"     , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(meta_rna    , "pipeline/1.1.metadata/meta_rna.txt"     , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(meta_subject, "pipeline/1.1.metadata/meta_subject.txt" , sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(meta_wgs    , "pipeline/1.1.metadata/meta_wgs.txt"     , sep = "\t", col.names = TRUE, row.names = FALSE)

# Write scRNA-seq metadata
meta_scrna           = fread("input/metadata/ipscore.scrna_metadata.txt", sep = "\t", header = TRUE, data.table = FALSE)
meta_scrna           = unique(meta_scrna[, c("biosample_source_id", "lab_frazer_udid")])
colnames(meta_scrna) = c("scrna_id", "udid")
meta_scrna           = meta_scrna[grepl(",", meta_scrna$scrna_id) == FALSE,]
meta_scrna           = merge(meta_scrna, meta_rna_ipscore[,c("udid", "subject_id", "subject_name", "assay_id")])

fwrite(meta_scrna, "pipeline/1.1.metadata/meta_scrna.txt" , sep = "\t", col.names = TRUE, row.names = FALSE)


