#build

library(here)

source(here("base", "src.R"))

#process data

source(here("process", "relative_abundance.R"))
source(here("process", "competencyCode.R"))
source(here("process", "weighted_prev.R"))
source(here("process", "data_format.R"))

#produce figures
source(here("scripts", "figures", "fig1.R"))
source(here("scripts", "figures", "fig2.R"))
source(here("scripts", "figures", "fig2_inset.R"))
source(here("scripts", "figures", "supp1.R"))
source(here("scripts", "figures", "supp2.R"))
source(here("scripts", "figures", "supp3_4.R"))
