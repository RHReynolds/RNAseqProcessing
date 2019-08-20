# RNAseq processing

## Aim
1. [Comparison of RNA QC run using either trimmomatic or fastp](comparison_trimmomatic_fastp/Comparison.md).
2. Provision of scripts for RNA-seq processing steps, from QC to alignment and quantification.

## Using the package

### Installing the package
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/RNAseqProcessing")
```

### Calling tools from the command line
The executables for all tools have already been downloaded to `/tools/`. To run all scripts below assumes that you are able to call various tools from the command line without having to point to their exact location on the server. To do this edit your .profile in your home directory:

```{bash, echo = T, eval = F}
cd 
nano .profile
```

Then ensure that each tool (fastp, STAR, etc.) has an `export PATH` line, which will tell bash to look in `/tools/` for commands. E.g. For fastp add: `export PATH="/tools/fastp/:$PATH"`. 

## Scripts for RNA-seq processing

 Script | Processing Step | Description | Author(s)
 ------ | --------------- | ----------- | ---------
 [prealignmentQC_fastp_PEadapters.R](QC/prealignmentQC_fastp_PEadapters.R) | Pre-alignment QC | This will perform fastp trimming, with adapter sequence auto-detection for PE data enabled, followed by fastQC and MultiQC. If you wish to specify adapters, this flag needs to be enabled. Script not yet produced. | DZ, KD & RHR
 [prealignmentQC_fastp_notrimming.R](QC/prealignmentQC_fastp_notrimming.R) | Pre-alignment QC | This will run fastp, but with trimming disabled, followed by fastQC and MultiQC. | DZ, KD & RHR
 [STAR_alignment_withReadGroups_multi2pass.R](alignment/STAR_alignment_withReadGroups_1pass.R) | Alignment | Performs STAR alignment, with the option of adding read groups if needed (this is important if you're planning to use you bams for later de-duplication with UMIs). By default, this script will perform 1st pass mapping. If users wish to use it for 2nd pass mapping, together with a file of filtered junctions, call the `--sj_file` flag. For details of alignment process, read the [alignment workflow](alignment/alignment.md). | DZ & RHR
  [STAR_splice_junction_merge.R](alignment/STAR_splice_junction_merge.R) | Alignment | Performs merging of SJ.out.tab files from 1st pass mapping, removes duplicated splice junctions (as determined by genomic location) and outputs one SJ.out.tab file with the genomic coordinates. Also has optional flag for filtering junctions by the number of samples they are present in. For details of alignment process, read the [alignment workflow](alignment/alignment.md). | RHR 
 [post_alignment_QC_RSeQC.R](QC/post_alignment_QC_RSeQC.R) | Post-alignment QC | Performs (i) sorting and indexing of .bam files using samtools and (ii) runs post-alignment QC, using RSeQC. For details, read the [alignment workflow](alignment/alignment.md). | DZ, KD & RHR
 [quantification_Salmon.R](quantification/quantification_Salmon.R) | Quantification | Performs mapping-based quantification of transcripts and genes (the latter is only if a transcript-to-gene map is provided). This script can be used following trimming, as it does not require aligned files. Instead, Salmon will perform quasi-mapping prior to quantification. The benefit of using Salmon for quantification is its speed and ability to correct for sequence-specific biases, GC-biases and positional biases. This script is adapted for paired-end reads. For more details, read the [quantification workflow](quantification/quantification.md). | RHR

## Example workflow
![](https://www.lucidchart.com/publicSegments/view/1973e61f-9de3-43eb-9261-d325df8e174c/image.png)
