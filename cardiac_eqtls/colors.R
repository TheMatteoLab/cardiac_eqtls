suppressPackageStartupMessages(library(colorspace    ))
suppressPackageStartupMessages(library(gameofthrones ))


study2color           = data.frame(study = c("ipscore", "gtex"), name = c("iPSCORE", "GTEx"), color = c("#0066CC", "#A52A56"), order = 2:1)
study2color           = study2color[order(study2color$order),]
rownames(study2color) = study2color$study
        
tissue2color           = data.frame(tissue    = c("arteria_aorta", "arteria_coronary", "ipsc_cvpc", "heart_atrium", "heart_ventricle"), 
                                    body_site = c("Aorta", "Coronary", "iPSC-CVPC", "Atrium", "Ventricle"), 
                                    color     = c("#8B636C", "#FFC0CB", "#0066CC", "#FF34B3", "#8B1C62"))

tissue2color[tissue2color$tissue == "ipsc_cvpc", "color"] = "#0066CC"
rownames(tissue2color) = tissue2color$tissue
tissue2color$order     = c(4,5,1,2,3)
tissue2color           = tissue2color[order(tissue2color$order),]

tissue2color3           = data.frame(tissue   = c("arteria", "ipsc_cvpc", "heart"  , "adult"  ), 
                                    body_site = c("Arteria", "iPSC-CVPC", "Heart"  , "Adult"  ), 
									color     = c("#2E8B57", "#0066CC"  , "#A52A56", "#A52A56"))

#tissue2color3[tissue2color3$tissue == "ipsc_cvpc", "color"] = "#FA8072"
rownames(tissue2color3) = tissue2color3$tissue
