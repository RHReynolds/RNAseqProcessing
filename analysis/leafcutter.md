# Table of contents
1. [Background](#background)
2. [Running the code](#running_the_code)
    1. [Example code](#example_code)
        1. [Differential splicing with multiple pairwise comparisons](#diff_splicing_multi)
        2. [Leafviz with multiple pairwise comparisons](#leafviz_multi)
        3. [Running Leafviz from inside RStudio](#leafviz_RStudio)
3. [Results](#Results)

# Background <a name="background"></a>

- Leafcutter reference: https://www.nature.com/articles/s41588-017-0004-9
- `leafcutter` quantifies intron excision ratios building a splicing graph with all split reads that have a shared donor or acceptor site sitting within one cluster
- `leafcutter` operates in 3 steps: 
    1. Converting aligned .bam files to .junc format detailing all the junctions of the each sample
    2. Generate a intron excision clusters and calculate the intron excision ratios 
    3. Perform differential splicing analysis 

```{r leafcutter fig 1a, echo=FALSE}
knitr::include_graphics("/home/rreynolds/projects/Aim2_PDsequencing_wd/raw_data/misc/Screenshot 2019-03-07 at 16.16.28.png")
```

- There are several filtering steps and parameters that are worth bearing in mind when using `leafcutter`. These include:
    - To identify clusters of introns, split reads that map with **>= 6 nt** into each exon are extract from aligned.bam files.
    - For each intron clusser, LeafCutter iteratively removes introns supported by 
        1. < number of reads across all samples (default: 30)
        2. OR < proportion of reads (default: 0.1%) of total number of intronic read counts for entire cluster
    - When performing differential splicing analyses, tests are only performed for:
        1. Introns detected (i.e. have >=1 corresponding spliced read) in >= select number of samples (default: 5)
        2. Clusters found in each group are detected in >= select number of individuals (default: 3), with 20 spliced reads supporting introns in the cluster.
    - These filters are customisable as optional parameters.
- **NOTE:** The documentation below provides an overview of the steps involved in Leafcutter differential splicing analyses, but not a step-by-step guide with commands. For a detailed tutorial on how to run Leafcutter, please refer to: http://davidaknowles.github.io/leafcutter/. The documentation below also highlights where `RNAseqProcessing` functions/scripts written for the Leafcutter pipeline can be used (these sections are preceded by the label, **`RNAseqProcessing`**), and where Leafcutter reference files can be found on the server.


# Running the code <a name="running_the_code"></a>

- The Leafcutter package can be found in the following directory on the server: `/tools/leafcutter/`. 
- Steps include:
    1. Generating .junc input files.
        1. By default, a .bam file is required for this step.
        2. **RNAseqProcessing:** If you have run STAR multi2pass alignment, you will have SJ.out.tab files, with the necessary information for Step 2. Instead of using .bam files, SJ.out.tab files can be formatted using the [`convert_STAR_SJ_to_junc()`](..R/leafcutter_functions.R) function, which can be called in RStudio.
    2. Clustering introns, using the `leafcutter_cluster.py` script in the Leafcutter package.
    3. Differential splicing analyses.
        1. Depending on the ensembl version used, this may require generation of a Leafcutter-appropriate exon file. This can be generated from a .gtf file using the `gtf_to_exons.R`script in the Leafcutter package. Before generating your own, check the following directory `/data/references/ensembl/gtf_gff3/` to see if a `leafcutter` directory with the necessary files already exists in the ensembl version required.
        2. **`RNAseqProcessing`** LeafCutter's differential splicing analyses currently only support pairwise comparisons. For each pairwise comparison it requires a group file to specify which samples belong to which group. Thus, if a grouping variable contains > 2 groups, multiple pairwise comparisons must be made and multiple group files generated. The [`create_group_files_multi_pairwisecomp()`](..R/leafcutter_functions.R) function can be used to do this. It will identify the comparisons, based on an inputted grouping column, and will output separate group .txt files for each group comparison combination.
        3. **`RNAseqProcessing`** With multiple comparisons, the differential splicing script provided by Leafcutter (`leafcutter_ds.R`) will have to be looped over the multiple comparisons. With the [leafcutter_ds_multi_pairwise.R](leafcutter_ds_multi_pairwise.R) script, which serves a wraparound for the original leafcutter script, this is possible.
    4. Visualise with Leafvis.
        1. For visualisation, Leafcutter requires a number of files, which can be produced from a .gtf file, using their provided `gtf2leafcutter.pl` script. Before generating your own, check the following directory `/data/references/ensembl/gtf_gff3/` to see if a `leafcutter` directory with the necessary files already exists in the ensembl version required.
        2. Results need to be prepared for use with Leafviz, which is done using the `prepare_results.R` script provided by Leafcutter.
        3. **`RNAseqProcessing`** To format the results of multiple pairwise comparisons requires looping across the various pairwise comparisons and running the `prepare_results.R` for each individual pairwise comparison. This is what the [`leafviz_multi_pairwise.R`](leafviz_multi_pairwise.R) script does.
            - Note: Using the `leafviz_multi_pairwise.R` script assumes use of the `leafcutter_ds_multi_pairwise.R` script and the `create_group_files_multi_pairwisecomp()` function. This ensures that all files needed are named consistently. That is, the (i) `_cluster_signficance.txt`, (ii) `_effect_sizes.txt` and (iii) `_group_file.txt` for an individual pairwise comparison all have the same prefix.  
        4. **`RNAseqProcessing`** Leafviz data can be viewed directly from RStudio, using the [`run_leafviz()`](..R/leafcutter_functions.R) function.

## Example code <a name="example_code"></a>

### Differential splicing with multiple pairwise comparisons <a name="diff_splicing_multi"></a>
```{bash run differential splicing leafcutter , echo = T, eval = F}
nohup Rscript /home/rreynolds/packages/RNAseqProcessing/analysis/leafcutter_ds_multi_pairwise.R \
/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/intron_clustering/tissue_polyA_test_diseasegroups_perind_numers.counts.gz \
/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/group_files/ \
--output_prefix=/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/ \
--max_cluster_size=Inf \
--min_samples_per_intron=5 \
--min_samples_per_group=3 \
--min_coverage=20 \
--timeout=30 \
--num_threads=15 \
--exon_file=/data/references/ensembl/gtf_gff3/v97/leafcutter/Homo_sapiens.GRCh38.97_LC_exon_file.txt.gz \
&>/home/rreynolds/projects/Aim2_PDsequencing_wd/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_leafcutter_ds.log&

```

### Leafviz with multiple pairwise comparisons <a name="leafviz_multi"></a>
```{bash run prepare_results , echo = T, eval = F}
nohup Rscript /home/rreynolds/packages/RNAseqProcessing/analysis/leafviz_multi_pairwise.R \
/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/intron_clustering/tissue_polyA_test_diseasegroups_perind_numers.counts.gz \
/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/ \
/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/ \
/data/references/ensembl/gtf_gff3/v97/leafcutter/Homo_sapiens.GRCh38.97 \
--output_dir=/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/leafviz/ \
--group_file_dir=/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/diff_splicing/group_files/ \
--FDR=0.05 \
&>/home/rreynolds/projects/Aim2_PDsequencing_wd/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_leafviz_prepare_results.log&

```

### Running Leafviz from inside RStudio <a name="leafviz_RStudio"></a>
```{r run Leafviz, echo = T, eval=F}

RNAseqProcessing::run_leafviz(leafviz_dir = "/home/rreynolds/packages/leafcutter/leafviz/", 
                              results_filename = "/home/rreynolds/projects/Aim2_PDsequencing_wd/results/leafcutter/leafviz/Control_vs_DLB.Rda")

```

# Results <a name="Results"></a>
If multiple pairwise comparisons have been performed, it is important that the p-values are corrected appropriately. This can be done by reading in all `_cluster_significance.txt` files and applying R's `p.adjust()` to all tests that were successfully tested (i.e. `status == "Success"`). 
