# Table of contents
1. [Download + installation](#download)
2. [Generating genome indexes](#genome_index)
    1. [Downloading the .fasta file](#downloading_fasta)
    2. [Downloading the .gtf file](#downloading_gtf)
    3. [Generating the genome index](#generate_genome_index)
3. [Mapping FASTQs with STAR](#mapping_w_star)
4. [Multi-sample 2-pass mapping](#multi_2pass)
5. [Post-alignment QC](#post_align_QC)
    1. [Sorting and indexing .bam files](#sort_and_index)
    2. [Converting gtf to bed format](#gtf_to_bed)
    3. [RSeQC](#RSeQC)
    4. [Running the post-alignment QC script](#running_post_align_QC)
6. [Multi-QC](#multiQC)

# Download + installation <a name="download"></a>

- Star downloaded from https://github.com/alexdobin/STAR/releases. Current version is 2.7.0a.
- Placed and unzipped into "/tools/STAR/". These come pre-compiled and the binary executable files are found in "/bin".
- Add the binary executable to PATH.

```{bash download/install STAR, echo = T, tidy = T, eval = F}
wget https://github.com/alexdobin/STAR/archive/2.7.0a.zip
unzip 2.7.0a.zip
export PATH=$PATH:/tools/STAR/STAR-2.7.0a/bin/Linux_x86_64/
```

# Generating genome indexes <a name="genome_index"></a>

- This generates a genome index that is used by STAR for mapping reads, improving the efficiency of alignment by not having to search across the entire genome, instead grouping regions of similar sequence into indexes.
- Requires as an input:
    1. The .fasta file detailing the genome
    2. A .gtf file detailing the annotations (check whether annotation has since been updated)
- **Note:** As users might have varying read lengths they may have to generate a new genome index specific to their read length.
    
## Downloading the .fasta file <a name="downloading_fasta"></a>
- The .fasta for hg38, ensembl v97 was downloaded using below command in the directory `/data/references/fasta`.
- In general, as long as the build does not change this file can be used repeatedly to generate genome indexes.

```{bash Download ensembl hg38 DNA reference, echo = T, tidy = F, eval = F}
wget ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
- Version number added to file name to distinguish from current versions in folder.

## Downloading the .gtf file <a name="downloading_gtf"></a>
- The latest version can be downloaded from: https://www.ensembl.org/info/data/ftp/index.html. 
- We have been storing various versions in the directory: `/data/references/ensembl/gtf_gff3/`. Choose the version suitable to your purposes.
- As an example, version 97 was downladed using the following command in the directory `/data/references/ensembl/gtf_gff3/v97`:

```{bash Download gtf, echo = T, tidy = F, eval = F}
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
gunzip Homo_sapiens.GRCh38.97.gtf.gz
```

## Generating the genome index <a name="generate_genome_index"></a>
- Below are the parameters used to generate the genome index: 
    + *runMode* - tells star we are generating an index and not mapping reads
    + *genomeDir* - output directory
    + *genomeFastaFiles* - fasta of genome
    + *sjdbGTFfile* - gtf file describing gene annotation 
    + *sjdbOverhang* - speciﬁes the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equalto the ReadLength-1, where ReadLength is the length of the reads. *In other words, this should be optimised to your read length.*
- Some genome indices have already been generated for different versions of ensembl. Please see the directory `/data/STAR_data/`. Within each version folder, you will find additional folders named `sjdbOverhang_.*`, with `.*` specifying the overhang used. 
- If you need to generate a new one, example code is shown below.

```{bash STAR genome index generation, echo = T, tidy = T, eval = F}
STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/STAR_data/genome_index_hg38_ens_v97/sjdbOverhang_99 \
--genomeFastaFiles /data/references/fasta/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa \
--sjdbGTFfile /data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf \
--sjdbOverhang 99
```

# Mapping FASTQs with STAR <a name="mapping_w_star"></a>

- For a sample of ~100mill 100bp paired end reads, using 15 threads, star takes about 15 minutes to finish aligning
- There are several parameters that can be set in STAR. Within both STAR alignment scripts, some of these arguments can be modified when calling the R script from the command line. Others are set, and have been set within the `get_star_parameters_set()` function in the [alignment_CommonFunctions.R](../R/alignment_CommonFunctions.R) script.  
- Modifiable arguments that are provided in the command line call, include:
```{bash, echo = T, tidy = T, eval = F}
  star_cmd <- str_c("STAR",
                    " --runThreadN ", threads_STAR,
                    " --genomeDir ", genome_index_path,
                    " --sj_file ", sj_path,
                    " --readFilesIn ", fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_paired[2],
                    " --readFilesCommand zcat ", # because fastq's are zipped
                    "--outFileNamePrefix ",  str_c(output_path, "/", sample_name_to_filter, "_ "))
```
- Set parameters include:
```{bash, echo = T, tidy = T, eval = F}

get_star_parameters_set <-function(){

  return(str_c(" --outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
               "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM,
               "--outFilterType BySJout ", # removes spurious split reads
               "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
               "--outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
               "--outFilterMismatchNoverReadLmax 0.04 ", # max number of mismatches per pair relative to read length. As per current ENCODE options.
               "--alignIntronMin 20 ", # min intron length. As per ENCODE options.
               "--alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
               "--alignMatesGapMax 1000000 ", # max gap between pair mates. As per ENCODE options.
               "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
               "--alignSJDBoverhangMin 3 " # minimum annotated split read anchor. Default is 3.
  ))

}

```
- **If set parameters require changing, please create a pull request with argumentation included as to why the lab should be changing these parameters.**
- STAR default parameters: https://github.com/alexdobin/STAR/blob/master/source/parametersDefault 
- STAR manual with details of more parameters: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- Many of the parameters set below were set according to the ENCODE standard options for long RNA-seq pipeline (as per [STAR manual 2.7.1a](STARmanual.pdf))
- To perform mapping use [STAR_alignment_withReadGroups_multi2pass.R](STAR_alignment_withReadGroups_multi2pass.R). Call the script using: `Rscript /path/to/script/STAR_alignment_withReadGroups_multi2pass.R -h`. The `-h` flag will list the required inputs and optional arguments.
    - It is important that required inputs are supplied in the correct order (i.e. fastq_dir_paths, genome_index_path, output_path).
    - Optional arguments can be called using flags. To see options, call the script with `-h` flag.
    - As an example, this script was run using the following arguments.
    
```{bash, echo = T, tidy = T, eval = F}
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/alignment/STAR_alignment_withReadGroups_multi2pass.R \
/data/RNAseq_PD/tissue_polyA_samples/QC/fastp \
/data/STAR_data/genome_index_hg38_ens_v97/sjdbOverhang_99 \
/data/RNAseq_PD/tissue_polyA_samples/STAR \
--sample_prefix=NM...._ \
--sample_suffix=_S.* \
--read_groups=/home/rreynolds/projects/Aim2_PDsequencing/data/Flowcell_info.txt \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_STAR_1pass.log&
```

# Multi-sample 2-pass mapping <a name="multi_2pass"></a>
- For sensitive novel junction discovery, 2-pass mapping is recommended. This can be run using:
    1. Multi-sample 2-pass mapping: implemented in this package.
    2. Per-sample 2-pass mapping: currently not implemented in this package.
- Multi-sample 2-pass mapping, is performed in three steps:
    1. Run 1st mapping pass for all samples with "usual" parameters using [STAR_alignment_withReadGroups_multi2pass.R](STAR_alignment_withReadGroups_multi2pass.R).
    2. Merge all SJ.out.tab files from all samples and remove duplicated junctions using [STAR_splice_junction_merge.R](STAR_splice_junction_merge.R). Note that this script also has an optional flag ```-f```, which adds a filter for the number of samples a junction must be found in to be included in the final merged file. E.g. if set to 5, then only junctions found in 5 or more samples will be included.
    3. Run 2nd mapping pass for all samples with "usual" parameters, and the addition of the merged SJ.out.tab file in the `--sj_file` flag. This is implemented in [STAR_alignment_withReadGroups_multi2pass.R](STAR_alignment_withReadGroups_multi2pass.R).
    
- As an example, step 2 and 3 were performed as a continuation of the commands in the last section, using the following commands:

```{bash, echo = T, tidy = T, eval = F}
# Merge SJ.out.tab files
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/alignment/STAR_splice_junction_merge.R \
/data/RNAseq_PD/tissue_polyA_samples/STAR/SJ_out_1pass \
-o /data/RNAseq_PD/tissue_polyA_samples/STAR/ \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_SJ_filtering.log&

# Run 2nd pass mapping
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/alignment/STAR_alignment_withReadGroups_multi2pass.R \
/data/RNAseq_PD/tissue_polyA_samples/QC/fastp \ 
/data/STAR_data/genome_index_hg38_ens_v97/sjdbOverhang_99 \
/data/RNAseq_PD/tissue_polyA_samples/STAR \
--sample_prefix=NM...._ \
--sample_suffix=_S.* \
--read_groups=/home/rreynolds/projects/Aim2_PDsequencing/data/Flowcell_info.txt \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_STAR_2pass.log&
--sj_file=/data/RNAseq_PD/tissue_polyA_samples/STAR/merged_junctions.SJ.out.tab \

```
# Post-alignment QC <a name="post_align_QC"></a>

## Sorting and indexing .bam files <a name="sort_and_index"></a>
- Many downstream applications require sorting and indexing of the .bam files. While STAR should, in theory, sort by co-ordinate, files did not seem to work with RSeQC, therefore built into [post_alignment_QC_RSeQC.R](../QC/post_alignment_QC_RSeQC.R) script a sorting and indexing using `samtools`. 
- **Note:** once the script has sorted and indexed the original .bam file outputted by STAR, the original STAR .bam file will be removed, such that only the samtools-sorted .bam file remains. 

## Converting gtf to bed format <a name="gtf_to_bed"></a>
- Several of the RSeQC modules require a reference gene model in the bed format. 
- Several already available on the server. Please check following directory: `/data/references/ensembl/bed/`
- If unavailable, can convert a .gtf to .bed format using the following commands:

```{bash, eval = F, echo = T}
/tools/ea-utils/ExpressionAnalysis-ea-utils-bd148d4/clipper/gtf2bed \
/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf \
> /data/references/ensembl/bed/v97/ensembl_GRCh38_v97.bed
```
## RSeQC <a name="RSeQC"></a>
- Used RSeQC for QC following alignment. 
- RSeQC Modules implemented include:
    - *geneBody_coverage.py*: calculates the RNA-seq reads coverage over gene body.
    - *clipping_profile.py*: estimates clipping profile of RNA-seq reads from .bam or .sam file. 
    - *inner_distance.py*: calculates the inner distance between read pairs.
    - *junction_annotation.py*: this compares detected splice junctions to a reference gene model.
    - *junction_saturation.py*: this module checks for saturation by resampling 5%, 10%, 15%, ..., 95% of total alignments from the .bam file and then detects splice junctions from each subset and compares them to the reference gene model.
    - *mismatch_profile.py*: calculates the distribution of mismatches across reads 
    - *read_distribution.py*: calculates how mapped reads were distributed over genomic feature e.g. exons, introns, intergenic, etc.
    - *read_duplication.py*
    - *read_GC.py*
    - *bam_stat.py*: summarises mapping statistics of .bam file
- For more details refer to: http://rseqc.sourceforge.net/

## Running the post-alignment QC script <a name="running_post_align_QC"></a>
As an example:

```{bash, echo = T, tidy = T, eval = F}
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/QC/post_alignment_QC_RSeQC.R \
/data/RNAseq_PD/tissue_polyA_samples/STAR \
/data/RNAseq_PD/tissue_polyA_samples/QC/ \
/data/references/ensembl/bed/v97/ensembl_GRCh38_v97.bed \
100 \
--sample_suffix=_Aligned.sortedByCoord.out.bam \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_nuclear_totalRNA_postalignment_QC.log&
```

# Multi-QC <a name="multiQC"></a>
Following on from post-alignment QC, it is possible to gather a full report of all QC, including pre-alignment QC. As an example:

```{bash, echo = T, tidy = T, eval = F}
multiqc /data/RNAseq_PD/tissue_polyA_samples/ \
-o /data/RNAseq_PD/tissue_polyA_samples/QC/multiqc/ \
--ignore /data/RNAseq_PD/tissue_polyA_samples/raw_data/ \
-n PD_tissue_polyA_full_report
```
