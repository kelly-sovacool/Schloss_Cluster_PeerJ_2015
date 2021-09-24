Why does VSEARCH perform worse than in Pat’s 2015 PeerJ paper?
================
2021-09-24

``` r
library(glue)
library(here)
library(tidyverse)
theme_set(theme_bw())
```

Pat had MCC’s of \~`0.7` for the MiSeq mouse dataset with VSEARCH DGC.
I’m getting MCC’s of `0.5781` in the OptiFit analysis repo with the same
dataset.

## My Snakefile

``` r
dat <- read_tsv(here('results/all_sensspec.tsv'))
dat %>% 
    ggplot(aes(dataset, mcc)) +
    geom_point() +
    facet_grid(mothur ~ vsearch) +
    ylim(0, 1)
```

![](/Users/kelly/projects/schloss-lab/Schloss_Cluster_PeerJ_2015/scratch/figures/dat-1.png)<!-- -->

## Pat’s Makefile

But when I run Pat’s makefile with mothur v1.37.0, get an MCC of `0.78`.
If I run it with mothur v.1.46.1, I get an MCC of `0`.
