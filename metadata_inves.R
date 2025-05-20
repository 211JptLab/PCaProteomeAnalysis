library(readr)
metadata <- read_delim("~/Downloads/metadata.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
library(dplyr)
table(metadata$Dataset)

