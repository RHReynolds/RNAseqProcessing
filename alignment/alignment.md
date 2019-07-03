### Download + installation 

- downloaded star from https://github.com/alexdobin/STAR/releases
- placed and unzipped into "/tools/STAR/"
- These come pre-compiled and the binary executable files are found in "/bin"
- Then the binary executable was added to PATH

```{r download/install STAR, echo = T, tidy = T, eval = F}

wget https://github.com/alexdobin/STAR/archive/2.7.0a.zip
unzip 2.7.0a.zip
export PATH=$PATH:/tools/STAR/STAR-2.7.0a/bin/Linux_x86_64/

```

<br> 

### Generating genome indexes

- This generates a genome index that is used by STAR for mapping reads, improving the efficiency of alignment by not having to search across the entire genome, instead grouping regions of similar sequence into indexes
- Requires as an input the .fasta file detailing the genome, a .gtf file detailing the annotations 
- The .fasta for hg38 was downloaded using command in the dir "/data/references/fasta": 

```{r Download ensembl hg38 reference, echo = T, tidy = F, eval = F}

rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-95/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

```

- below are the parameters used to generate the genome index: 
    + *runMode* - tells star we are generating an index and not mapping reads
    + *genomeDir* - output dir
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

<br>

### Mapping FASTQs with STAR

- Below is the annotated code used to map fastqs with STAR
- For a sample of ~100mill 100bp paired end reads, using 15 threads, star takes about 15 minutes to finish aligning
- Below parameters: 
    + *runThreadN* - number of threads to use, testing this mrserver doesn't allow you to go past 18
    + *genomeDir* - genome index directory path created above
    + *readFilesIn* - the two input fastq's (paired end) that you would like to align
    + *readF
- STAR manual with details of more parameters (http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf)

```{r STAR funning, echo = T, tidy = T, eval = F}

  system(command = str_c("/tools/STAR/STAR-2.7.0a/bin/Linux_x86_64/STAR", 
                         " --runThreadN ", threads_STAR, 
                         " --genomeDir ", genome_index_path, 
                         " --readFilesIn ", fastq_per_sample_paths_trimmed_paired[1], " ", fastq_per_sample_paths_trimmed_paired[2],  
                         " --readFilesCommand  zcat ", # because fastq's are zipped
                         "--outFilterType BySJout ", # removes spurious split reads
                         "--outFilterMultimapNmax 1 ", # only allows reads to be mapped to one position in the genome
                         "--alignSJoverhangMin 8 ", # minimum unannotated split read anchor 
                         "--alignSJDBoverhangMin 3 ", # maximum unannotated split read anchor 
                         "--outFilterMismatchNmax 3 ", # max num mismatches between pairs (was 2 for )
                         "--alignIntronMin 20 ", # min intron length (currently from ensembl its )
                         "--alignIntronMax 1000000 ", # max intron length (currently from ensembl its 1,097,903 from KCNIP4)
                         "--alignMatesGapMax 1000000 ", # max gap between pair mates
                         "--outSAMtype BAM SortedByCoordinate ", # output as a sorted BAM
                         "--outFileNamePrefix ",  str_c(output_path, "/", sample_name_to_filter, "_")))

```
