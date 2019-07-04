** WORK IN PROGRESS **

# RNAseq processing

## Aim
1. [Comparison of RNA QC run using either trimmomatic or fastp](comparison_trimmomatic_fastp/Comparison.md).
2. Provision of scripts for QC and alignment.

## Using the package

### Installing the package
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/RNAseqProcessing", auth_token = "")
```

As this is a private repository, you will have to generate a [personal access token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line), and insert this into the ```auth_token``` argument. **Remember to save this token, as you may need it to access other private repositories.**

### Calling tools from the command line
The executables for all tools have already been downloaded to `/tools/`. To run all scripts below assumes that you are able to call various tools from the command line without having to point to their exact location on the server. To do this edit your .profile in your home directory:

```{bash, echo = T, eval = F}
cd 
nano .profile
```

Then ensure that each tool (fastp, STAR, etc.) has an `export PATH` line, which will tell bash to look in `/tools/` for commands. E.g. For fastp add: `export PATH="/tools/fastp/:$PATH"`. 

## Scripts for QC and alignment

 Script | Processing Step | Description | Author(s)
 ------ | --------------- | ----------- | ---------
 [prealignmentQC_fastp_PEadapters.R](QC_RNAseq_samples/prealignmentQC_fastp_PEadapters.R) | Pre-alignment QC | This will perform fastp trimming, with adapter sequence auto-detection for PE data enabled. If you wish to specify adapters, this flag needs to be enabled. Script not yet produced. | DZ, KD & RHR
 [prealignmentQC_fastp_notrimming.R](QC_RNAseq_samples/prealignmentQC_fastp_notrimming.R) | Pre-alignment QC | This will run fastp, but with trimming disabled. | DZ, KD & RHR
 [STAR_alignment_1pass.R](alignment/STAR_alignment.R) | Alignment | Performs STAR alignment, using 1-pass mapping. This script can be combined with [STAR_alignment_multi2pass.R](STAR_alignment_multi2pass.R) for multi-sample 2-pass mapping. For details of alignment process, read the [alignment workflow](alignment/alignment.md). **N.B. This script is currently optimised for 100bp reads.** | DZ
 [STAR_alignment_withReadGroups_1pass.R](alignment/STAR_alignment_withReadGroups_multi2pass.R) | Alignment | Performs STAR alignment, using 1-pass mapping, and adds read groups, which can be used for later de-duplication with UMIs. This script can be combined with [STAR_alignment_multi2pass.R](STAR_alignment_multi2pass.R) for multi-sample 2-pass mapping. For details of alignment process, read the [alignment workflow](alignment/alignment.md). **N.B. This script is currently optimised for 100bp reads.** | DZ & RHR
 [STAR_alignment_multi2pass.R](alignment/STAR_alignment_multi2pass.R) | Alignment | Performs multi-sample 2-pass mapping. This is recommended for a study with multiple samples, and requires that 1-pass mapping has been performed, all junctions from each SJ.out.tab in the 1st pass have been collated, and that some filtering has occurred on these junctions (if no filtering is performed, the number of junctions will overload the server). Filtering will be study-dependent. | RHR 

