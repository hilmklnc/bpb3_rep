<!-- 
title: 'BPB3: BayesPI-BAR version 3 package User Guide'
author: "Junbai Wang, Mingyi Yang, Magnar Bjoras, ..."
date: "`r format(Sys.time(), '%d %B, %Y')`" 
output: 
 html_document: 
 df_print: paged 
 word_document: default
-->
---
title: 'BPB3: </br>
        Bayesian method for Protein–DNA Interaction with Binding Affinity Ranking version 3 (BayesPI-BAR3)'
---

<!--- R markdown set up --->
<!--
```{r setup, include=T, echo=FALSE, error=T, message=F}
knitr::opts_chunk$set(echo = TRUE)
library(rmarkdown)
library(stargazer)
library(reticulate)
reticulate::use_python("/Users/junbai/miniconda3/bin//python")
knitr::opts_chunk$set(engine.path = list(python = '/Users/junbai/miniconda3/bin//python'))
```
-->


# Install BPB3 package 
#
# Please send email to Dr. Junbai Wang (junbai.wang@medisin.uio.no) for requesting a free download of bpb3 pacakge. 
#
##  Please include your name, Institution/company name, Email of your institute or company (do not use gmail or yahoo email etc), and a brief description of the purpose for using bpb3 package.
##  Please use the subject title "bpb3 request "  for your email.
#
## Download and unzip package

```{bash,  echo=T}
# tar -zxf bpb3.tgz
# cd bpb3
# python setup.py install
```

## Tasks
```
Type "bpb3 --help" to list tasks in bpb3 package.
Then type "bpb3 <task> --help" to display help of each task.
````

```
bpb3 --help
usage:  bpb3 <task> [<args>]

     Tasks available for using:
         differential_expression        Predict differentialy expressed genes (DEG) based on 
                                        two group of samples.
         gene_regions           Extracts regions near transcription start sites of selected 
                                genes based on genCode gtf. 
         mussd                  Mutation filtering based on the Space and Sample Distribution 
                                - MuSSD. 
         highly_mutated_blocks  Find blocks with significantly more mutations than would be expected.
         bayespi_bar            BayesPI-BAR delta-dbA ranking computation for TF binding affinity 
                                affected by DNA mutation.
         choose_background_parameters   Selects parameters for mutation background computation.
         background_affinity_changes    Mutation background computation.
         affinity_change_significance_test      Significant test of TF binding affinity changes 
                                                between foreground and background affinity changes.
         parallel               Run commands from the given file in parallel.
         make_cluster4pwm       Make input PWM files for bpb3 based on clustered PWMs.
         bpb3selectedPWM        The second level analysis of bpb3 by using the top PWMs in TF ranking 
                                after the first level analysis of bpb3 based on the clustered PWMs. 
         run_pipeline           Run full bpb3 pipeline (e.g., the first level analysis of bpb3 if 
                                clustered PWMs are used in the calculation).
         clean_tmp              Clean temporary files from output folders.

     Tasks available for demo purpose:
         plot_result            Generate heatmaps for selected mutation blocks. (demo) 
         filter_results_by_gene_expression_cluster4pwm  Filter those TF whose expression is too 
                                                        low in clustered PWMs. (demo)
         filter_results_by_gene_expression              Filter those TFs whose expression is too 
                                                        low (e.g., RPKM<0.03). (demo)
         make_plots_cluster4pwm         Make heatmap plots for all significant mutation blocks 
                                        that affecting clustered PWMs. (demo)
         make_plots                     Make heatmap plots for all significant mutation blocks 
                                        that affecting PWMs. (demo)
         check_accuracy4cluster Check accuracy for 67 SNPs that based on clustered PWMs. (demo)
         check_accuracy         Check accuracy for 67 SNPs that based on original PWMs. (demo)
         filterDEG4bpb3         Filter bpb3 exported differential expression gene list by rratios. (demo)
         preprocess_icgc_data   Preprocess of ICGC data such as a folder contains files donor_*, 
                                specimen, simple_somatic*, exp_seq.tsv et al. (demo)
     
 BayesPI-BAR in Python3 - bpb3

 positional arguments:
  task        Pipeline task to run

 optional arguments:
  -h, --help  show this help message and exit
```

## Parallel options
```
Run bpb3 in a local computer with multiple processes: only needs four lines, "--user_cores", number of cores, "--max_nodes", and maximum number of nodes in "parallel_options.txt".

Run bpb3 in slurm:  add three new lines "--use_slurm" , "--slurm_account" , and slurm account name in "parallel_options.txt".
```

## BPB3 Demos and User Guide

[More information about demos in bpb3 package please refer to ->  https://bpb3.github.io/bpb3-web ](https://bpb3.github.io/bpb3-web/)

[bpb3 package can be downloaded from -> https://github.com/junbaiw/bpb3 ](https://github.com/junbaiw/bpb3)

## Readme in bpb3

*BayesPI-BAR3*  
Type "python setup.py install" to install bpb3  package on your local machine.        
file folders: final_demo, final_demo_data

## Reference papers

1. Junbai Wang, Kirill Batmanov, BayesPI-BAR: a new biophysical model for characterization of regulatory sequence variations, Nucleic Acids Research, Volume 43, Issue 21, 2 December 2015, Page e147. [Link of Paper 1](https://academic.oup.com/nar/article/43/21/e147/2468103)

2. Batmanov, K., Wang, W., Bjørås, M. et al. Integrative whole-genome sequence analysis reveals roles of regulatory mutations in BCL6 and BCL2 in follicular lymphoma. Sci Rep 7, 7040 (2017). [Link of Paper 2](https://www.nature.com/articles/s41598-017-07226-4)

3. Batmanov K, Delabie J, Wang J. BayesPI-BAR2: A New Python Package for Predicting Functional Non-coding Mutations in Cancer Patient Cohorts. Front Genet. 2019. [Link of Paper 3](https://www.frontiersin.org/articles/10.3389/fgene.2019.00282/full)

4. Farooq A, Trøen G, Delabie J, Wang J. Integrating whole genome sequencing, methylation, gene expression, topological associated domain information in regulatory mutation prediction: A study of follicular lymphoma. Comput Struct Biotechnol J. 2022 Mar 23;20:1726-1742. [Link of Paper 4](https://www.sciencedirect.com/science/article/pii/S2001037022000976)   

5. Mingyi Yang, Omer Ali, Magnar Bjørås, Junbai Wang. Identifying Functional Regulatory Mutation Blocks by Integrating Genome Sequencing and Transcriptome Data. iScience 2023 July 2. [link of Paper 5](https://doi.org/10.1016/j.isci.2023.107266)


## Notes

1. Demo1 and demo2 were used in papers 1 and 3.
2. Demo3, demo4, demo5, and demo6 were used in papers 2 and 3.
3. Demo7 was similar to the calculation in paper 4, but the 14 FL samples were replaced by an independent FL cohort from papers 2 and 3.
4. Under final_demo_data, the pwm folder contains 1772 human pwms that were used in papers 1, 2, 3, and 4.
5. Applications of clustered PWMs in demos 2, 4, and 5 are described in current paper 5.  

# End


