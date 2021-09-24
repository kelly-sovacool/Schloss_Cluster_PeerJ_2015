library(tidyverse)

read_tsv(snakemake@input[['tsv']]) %>%
    mutate(dataset = snakemake@wildcards[['dataset']],
           method = snakemake@wildcards[['method']],
           vsearch = snakemake@wildcards[['vver']],
           mothur = snakemake@wildcards[['mver']]
       ) %>%
    write_tsv(snakemake@output[['tsv']])