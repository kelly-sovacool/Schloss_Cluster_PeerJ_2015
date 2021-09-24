library(tidyverse)

snakemake@input[['tsv']] %>%
  map_dfr(read_tsv) %>%
  write_tsv(snakemake@output[['tsv']])
