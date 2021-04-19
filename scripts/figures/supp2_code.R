#Daniel Suh and Andrew Park
#6/3/20

#Figure 2 in rv_cc manuscript

library(here)
library(metacom)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))

community_mat %<>% column_to_rownames(., var = "siteID")

Imagine(t(community_mat), 
        fill = FALSE, 
        order = TRUE, 
        xlab = "", 
        ylab = "",
        xline = 5,
        yline = 4,
        sitenames = ,
        speciesnames = ,
        binary=T)

ordered_matrix <- OrderMatrix(community_mat)
rownames(ordered_matrix)

results <- Metacommunity(community_mat, binary = T)

#wasn't able to save output from Imagine as a plot and had to save using dialog within RStudio