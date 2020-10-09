#Script Originator: Daniel Suh

#This script creates figure 7


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)

data <- read_csv(here("data/weighted_prev_competence.csv"))
