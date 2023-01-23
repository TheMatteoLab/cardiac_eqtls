suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(gplots    ))
suppressPackageStartupMessages(library(plyr      ))
suppressPackageStartupMessages(library(dplyr     ))
suppressPackageStartupMessages(library(seqminer  ))

add_rownames = function(x) # add rownames to fread
{
	rownames(x) = x[,1]
	x[,1]       = NULL
	return(x)
}

values2color = function(x, colorscale = colorpanel(n = 100, low = "#ffeeee", high = cm.color), minval = NA, maxval = NA) # convert numeric to color scale
{
    if (is.na(minval) == TRUE){minval = min(x)}
    if (is.na(maxval) == TRUE){maxval = min(x)}
    if (length(x[is.na(x) == TRUE]) > 0){x[is.na(x) == TRUE] = minval}
    if (length(x[x < minval      ]) > 0){x[x < minval      ] = minval}
    if (length(x[x > maxval      ]) > 0){x[x > maxval      ] = maxval}

    y = round((x - minval) * length(colorscale) / (maxval - minval))
    if (length(y[y < 1]) > 0){y[y < 1] = 1}
    
    z = colorscale[y]
    
    return(z)
}


