** WORK IN PROGRESS **

## Download + installation 

- Star downloaded from https://github.com/alexdobin/STAR/releases. Current version is 2.7.0a.
- Placed and unzipped into "/tools/STAR/". These come pre-compiled and the binary executable files are found in "/bin".
- Add the binary executable to PATH.

```{r download/install STAR, echo = T, tidy = T, eval = F}

wget https://github.com/alexdobin/STAR/archive/2.7.0a.zip
unzip 2.7.0a.zip
export PATH=$PATH:/tools/STAR/STAR-2.7.0a/bin/Linux_x86_64/

```

## Generating genome indexes

- This generates a genome index that is used by STAR for mapping reads, improving the efficiency of alignment by not having to search across the entire genome, instead grouping regions of similar sequence into indexes.
- Requires as an input the .fasta file detailing the genome, a .gtf file detailing the annotations 
- The .fasta for hg38 was downloaded using command in the dir "/data/references/fasta": 

```{r Download ensembl hg38 reference, echo = T, tidy = F, eval = F}

rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-95/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

```

- Below are the parameters used to generate the genome index: 
    + *runMode* - tells star we are generating an index and not mapping reads
    + *genomeDir* - output directory
    + *genomeFastaFiles* - fasta of genome
    + *sjdbGTFfile* - gtf file describing gene annotation 
    + *sjdbOverhang* - STAR says this should be optimally (read length - 1) - can be seen as the maximum split read anchor length
    
```{r STAR genome index generation, echo = T, tidy = T, eval = F}

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/STAR_data/genome_index_hg38 \
--genomeFastaFiles /data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /data/references/ensembl/gtf_gff3/v95/Homo_sapiens.GRCh38.95.gtf \
--sjdbOverhang 149

```

## Mapping FASTQs with STAR

- For a sample of ~100mill 100bp paired end reads, using 15 threads, star takes about 15 minutes to finish aligning
- There are several parameters that can be set in STAR. Below is what has been used in [STAR_alignment_1pass.R](Star_alignment.R).
- STAR default parameters: https://github.com/alexdobin/STAR/blob/master/source/parametersDefault 
- STAR manual with details of more parameters: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

```{r STAR function, echo = T, tidy = T, eval = F}

  system(command = str_c("STAR",
                         " --runThreadN ", threads_STAR,
                         " --genomeDir ", genome_index_path,
                         " --readFilesIn ", fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_paired[2],
                         " --readFilesCommand  zcat ", # because fastq's are zipped
                         "--outFilterType BySJout ", # removes spurious split reads
                         "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
                         "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor
                         "--alignSJDBoverhangMin 3 ", # minimum annotated split read anchor
                         "--outFilterMismatchNmax 2 ", # max num mismatches between pairs. For 100 bp read, allow 2. For 150 bp reads, use 3.
                         "--alignIntronMin 20 ", # min intron length
                         "--alignIntronMax 1000000 ", # max intron length (currently from ensembl its 1,097,903 from KCNIP4)
                         "--alignMatesGapMax 1000000 ", # max gap between pair mates
                         "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM
                         "--outFileNamePrefix ",  str_c(output_path, "/", sample_name_to_filter, "_")))

```

## 2-pass mapping
- For sensitive novel junction discovery, 2-pass mapping is recommended. This can be run using:
    1. Multi-sample 2-pass mapping: this is currently implemented in [STAR_alignment_multi2pass.R](STAR_alignment_multi2pass.R)
    2. Per-sample 2-pass mapping: currently not implemented in a script.
