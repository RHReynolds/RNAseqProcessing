# RNAseq processing

## Aim
1. [Comparison of RNA QC run using either trimmomatic or fastp](Comparison_Trimmomatic_Fastp/Comparison.md).
2. Provision of scripts for QC and alignment.

## Using the package
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/RNAseqProcessing", auth_token = "")
```

As this is a private repository, you will have to generate a [personal access token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line), and insert this into the ```auth_token``` argument. **Remember to save this token, as you may need it to access other private repositories.**

## Scripts for QC and alignment

 Script | Processing Step | Description | Author(s)
 ------ | --------------- | ----------- | ---------
 [prealignmentQC_fastp_PEadapters.R](QC_RNAseq_samples/prealignmentQC_fastp_PEadapters.R) | Pre-alignment QC | This will perform fastp trimming, with adapter sequence auto-detection for PE data enabled. If you wish to specify adapters, this flag needs to be enabled. Script not yet produced. | DZ, KD & RHR
 [prealignmentQC_fastp_notrimming.R](QC_RNAseq_samples/prealignmentQC_fastp_notrimming.R) | Pre-alignment QC | This will run fastp, but with trimming disabled. | DZ, KD & RHR
