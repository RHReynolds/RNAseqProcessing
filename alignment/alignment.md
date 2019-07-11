** WORK IN PROGRESS **

## Download + installation 

- Star downloaded from https://github.com/alexdobin/STAR/releases. Current version is 2.7.0a.
- Placed and unzipped into "/tools/STAR/". These come pre-compiled and the binary executable files are found in "/bin".
- Add the binary executable to PATH.

```{bash download/install STAR, echo = T, tidy = T, eval = F}
wget https://github.com/alexdobin/STAR/archive/2.7.0a.zip
unzip 2.7.0a.zip
export PATH=$PATH:/tools/STAR/STAR-2.7.0a/bin/Linux_x86_64/
```

## Generating genome indexes

- This generates a genome index that is used by STAR for mapping reads, improving the efficiency of alignment by not having to search across the entire genome, instead grouping regions of similar sequence into indexes.
- Requires as an input:
    1. The .fasta file detailing the genome
    2. A .gtf file detailing the annotations (check whether annotation has since been updated)
- *Note:* As users might have varying read lengths they may have to generate a new genome index specific to their read length.
    
### Downloading the .fasta file
- The .fasta for hg38, ensembl v97 was downloaded using below command in the directory `/data/references/fasta`.
- In general, as long as the build does not change this file can be used repeatedly to generate genome indexes.

```{bash Download ensembl hg38 DNA reference, echo = T, tidy = F, eval = F}
wget ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
- Version number added to file name to distinguish from current versions in folder.

### Downloading the .gtf file
- The latest version can be downloaded from: https://www.ensembl.org/info/data/ftp/index.html. 
- We have been storing various versions in the directory: `/data/references/ensembl/gtf_gff3/`. Choose the version suitable to your purposes.
- As an example, version 97 was downladed using the following command in the directory `/data/references/ensembl/gtf_gff3/v97`:

```{bash Download gtf, echo = T, tidy = F, eval = F}
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
gunzip Homo_sapiens.GRCh38.97.gtf.gz
```

### Generating the genome index
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

## Mapping FASTQs with STAR

- For a sample of ~100mill 100bp paired end reads, using 15 threads, star takes about 15 minutes to finish aligning
- There are several parameters that can be set in STAR. Below is what has been used in [STAR_alignment_withReadGroups_1pass.R](STAR_alignment_withReadGroups_1pass.R).

```{r STAR function, echo = T, tidy = T, eval = F}
  system(command = str_c("STAR",
                           " --runThreadN ", threads_STAR,
                           " --genomeDir ", genome_index_path,
                           " --readFilesIn ", fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_paired[2],
                           " --readFilesCommand  zcat ", # because fastq's are zipped
                           "--outFileNamePrefix ",  str_c(output_path, "/", sample_name_to_filter, "_ "),
                           "--outReadsUnmapped Fastx ", # output in separate fast/fastq files the unmapped/partially-mapped reads
                           "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM
                           "--outSAMattrRGline ", str_c("ID:", read_group$full_name, " PU:", read_group$PU, " SM:", read_group$sample_name, " PL:Illumina LB:xxx "), # SAM/BAM read group line
                           "--outFilterType BySJout ", # removes spurious split reads
                           "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
                           "--outFilterMismatchNmax 999 ", # Maximum number of mismatches per pair. Large numbers switch off filter. Instead we filter by "--outFilterMismatchNoverReadLmax".
                           "--outFilterMismatchNoverReadLmax 0.04 ", # max number of mismatches per pair relative to read length. As per current ENCODE options.
                           "--alignIntronMin 20 ", # min intron length. As per ENCODE options.
                           "--alignIntronMax 1000000 ", # max intron length. As per ENCODE options (currently from ensembl its 1,097,903 from KCNIP4).
                           "--alignMatesGapMax 1000000", # max gap between pair mates. As per ENCODE options.
                           "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor. As per ENCODE options.
                           "--alignSJDBoverhangMin 3" # minimum annotated split read anchor. Default is 3.
                           ))
```

- STAR default parameters: https://github.com/alexdobin/STAR/blob/master/source/parametersDefault 
- STAR manual with details of more parameters: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- Many of the parameters set below were set according to the ENCODE standard options for long RNA-seq pipeline (as per [STAR manual 2.7.1a](STARmanual.pdf))
- To perform mapping use [STAR_alignment_withReadGroups_1pass.R](STAR_alignment_withReadGroups_1pass.R). Call the script using: `Rscript /path/to/script/STAR_alignment_withReadGroups_1pass.R -h`. The `-h` flag will list the required inputs and optional arguments.
    - It is important that required inputs are supplied in the correct order (i.e. fastq_dir_paths, genome_index_path, output_path).
    - Optional arguments can be called using flags. To see options, call the script with `-h` flag.
    - As an example, this script was run using the following arguments.
    
```{bash, echo = T, tidy = T, eval = F}
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/alignment/STAR_alignment_withReadGroups_1pass.R \
/data/RNAseq_PD/tissue_polyA_samples/QC/fastp \
/data/STAR_data/genome_index_hg38_ens_v97/sjdbOverhang_99 \
/data/RNAseq_PD/tissue_polyA_samples/STAR \
--sample_prefix=NM...._ \
--sample_suffix=_S.* \
--read_groups=/home/rreynolds/projects/Aim2_PDsequencing/data/Flowcell_info.txt \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_STAR_1pass.log&
```

## Multi-sample 2-pass mapping
- For sensitive novel junction discovery, 2-pass mapping is recommended. This can be run using:
    1. Multi-sample 2-pass mapping: implemented in this package.
    2. Per-sample 2-pass mapping: currently not implemented in this package.
- Multi-sample 2-pass mapping, is performed in three steps:
    1. Run 1st mapping pass for all samples with "usual" parameters using [STAR_alignment_withReadGroups_1pass.R](STAR_alignment_withReadGroups_1pass.R).
    2. Merge all SJ.out.tab files from all samples and remove duplicated junctions using [STAR_splice_junction_merge.R](STAR_splice_junction_merge.R).
    3. Run 2nd mapping pass for all samples with "usual" parameters, and the addition of the merged SJ.out.tab file in the `--sjdbFileChrStartEnd` flag. This is implemented in [STAR_alignment_multi2pass.R](STAR_alignment_multi2pass.R).
    
- As an example, step 2 and 3 were performed as a continuation of the commands in the last section, using the following commands:

```{bash, echo = T, tidy = T, eval = F}
# Merge SJ.out.tab files
Rscript /home/rreynolds/projects/RNAseqProcessing/alignment/STAR_splice_junction_merge.R /data/RNAseq_PD/tissue_polyA_samples/STAR

# Run 2nd pass mapping
nohup Rscript \
/home/rreynolds/projects/RNAseqProcessing/alignment/STAR_alignment_multi2pass.R \
/data/RNAseq_PD/tissue_polyA_samples/QC/fastp \ 
/data/STAR_data/genome_index_hg38_ens_v97/sjdbOverhang_99 \
/data/RNAseq_PD/tissue_polyA_samples/STAR \
/data/RNAseq_PD/tissue_polyA_samples/STAR/all_samples_non_duplicated_junctions.SJ.out.tab \
--sample_prefix=NM...._ \
--sample_suffix=_S.* \
--read_groups=/home/rreynolds/projects/Aim2_PDsequencing/data/Flowcell_info.txt \
&>/home/rreynolds/projects/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_STAR_2pass.log&
```


