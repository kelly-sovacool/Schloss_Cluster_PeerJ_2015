***De novo* clustering methods out-perform reference-based methods for assigning 16S rRNA gene sequences to operational taxonomic units**
=======

**Background.** 16S rRNA gene sequences are routinely assigned to operational taxonomic units (OTUs) that are then used to analyze complex microbial communities. A number of methods have been employed to carry out the assignment of 16S rRNA gene sequences to OTUs leading to confusion over which method is the most rigorous. A recent study suggested that a clustering method should be selected based on its ability to generate stable OTU assignments that do not change as additional sequences are added to the dataset. In contrast, we contend that the ability of the method to properly represent the distances between the sequences is more important.

**Methods.** Our analysis implemented five *de novo* clustering algorithms including the single linkage, complete linkage, average linkage, abundance-based greedy clustering, distance-based greedy clustering, and the open and closed-reference methods. Using two previously published datasets we used the Matthew’s Correlation Coefficient (MCC) to assess the stability and quality of OTU assignments.

**Results.** The stability of OTU assignments did not reflect the quality of the assignments. Depending on the dataset being analyzed, the average linkage and the distance and abundance-based greedy clustering methods generated more robust OTUs than the open and closed-reference methods. We also demonstrated that for the greedy algorithms VSEARCH produced assignments that were comparable to those produced by USEARCH making VSEARCH a viable free and open source alternative to USEARCH. Further interrogation of the reference-based methods indicated that when USEARCH is used to identify the closest reference, the OTU assignments were sensitive to the order of the reference sequences because the reference sequences can be identical over the region being considered. More troubling was the observation that while both USEARCH and VSEARCH have a high level of sensitivity to detect reference sequences, the specificity of those matches was poor relative to the true best match.

**Discussion.** Our analysis calls into question the quality and stability of OTU assignments generated by the open and closed-reference methods as implemented in current version of QIIME. This study demonstrates that *de novo* methods are the most rigorous and that the quality of clustering assignments needs to be assessed for multiple methods to identify the optimal clustering method for a particular dataset.



Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw/         # raw data, will not be altered
    |  |- he/	       # direct replication of He et al. analysis
    |  |- schloss/     # how I would have processed Canadian soil data
    |  |- miseq/       # analysis of murine MiSeq data
    |  |- gg_13_8/     # analysis of QIIME reference database
    |  |- rand_ref/    # analysis of murine MiSeq data with randomized database
    |  +- process/     # cleaned data, will not be altered once created;
    |                  # will be committed to repo
    |
    |- code/           # any programmatic code
    |
    |- results         # all output from workflows and analyses
    |  |- figures/     # graphs, likely designated for manuscript figures
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
	|- papers/		   # files used to write, submit, and publish paper
	|  |- header.tex      # LaTeX header file for formatting paper
    |  |
	|  +- peerj_2015/    # original study
	|  |  |- Schloss_Cluster_PeerJ_2015.Rmd # executable Rmarkdown for this study
	|  |  |- Schloss_Cluster_PeerJ_2015.md  # Markdown (GitHub) version of the *Rmd file
	|  |  |- Schloss_Cluster_PeerJ_2015.pdf # PDF version of *.Rmd file
	|  |  |- peerj.csl       # CSL file for formatting bibliograph using PeerJ's format
	|  |  + references.bib  # bibtex formatted file of references
	|  |
	|  + msystems_2016 # response to Kopylova study
	|     |- Schloss_Cluster_mSystems_2016.Rmd # executable Rmarkdown for this study
	|     |- Schloss_Cluster_mSystems_2016.md  # Markdown (GitHub) version of the *Rmd file
	|     |- Schloss_Cluster_mSystems_2016.pdf # PDF version of *.Rmd file
	|     |- peerj.csl       # CSL file for formatting bibliograph using PeerJ's format
	|     + references.bib  # bibtex formatted file of references
	|
    |- Makefile        # executable Makefile for this study
    |
    +- LICENSE.md



Dependencies
------------
The following need to be installed and in the path...
* mothur (v.1.37.0)
* QIIME (v.1.9.1)
* UCLUST (v.6.1.544)
* VSEARCH (v.1.5.0)
* micca (OTUClust)
* R
    + jsonlite
    + knitr
    + Rcpp
    + plyr
    + ggplot2



These should be installed the specified folders
* vsearch (v.1.5.0 installed in code/vsearch)
* sumaclust (v.1.0.20 installed in code/sumaclust_v1.0.20)
* swarm (v.2.1.1 installed in code/swarm)
* NINJA-OPS (v.1.5.0 installed in code/NINJA-OPS
