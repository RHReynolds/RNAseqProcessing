# Table of contents
1. [Download + installation](#download)
2. [Running Salmon](#running_salmon)
3. [Creating a salmon index](#creating_salmon_index)
    a. [Example of creating a salmon index ](#example_salmon_index)
4. [Quantification](#quantification)
    a. [Determining the library type](#library_type)
    b. [Running Salmon quantification](#running_salmon_quant)
    c. [Loading the data post-quantification](#loading_post_quant_data)
        i. [Creating a named vector](#creating_named_vector)
        ii. [Creating a transcript-to-gene map](#transcript_to_gene_map)
    d. [Post-quantification QC checks](#QC_checks)

## Download + installation <a name="download"></a>

- Salmon downloaded from https://github.com/alexdobin/STAR/releases. Current version is 0.14.1.
- Placed and unzipped into `/tools/salmon/salmon-latest_linux_x86_64/`. These come pre-compiled and the binary executable files are found in `/bin`.
- Optional: add the binary executable to PATH.

```{bash download/install salmon, echo = T, tidy = T, eval = F}
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz
tar -xvzf salmon-0.14.1_linux_x86_64.tar.gz
export PATH=$PATH:/tools/salmon/salmon-latest_linux_x86_64/bin/
```

## Running Salmon <a name="running_salmon"></a>
- While Salmon can be run using previously aligned .bam files, it requires that these were aligned directly to the transcriptome rather than to the genome. There are ways of converting genome-aligned .bam files; however, given the speed of Salmon mapping, opted to perform quasi-mapping and quantification (otherwise known as mapping-based mode in Salmon documentation - https://salmon.readthedocs.io/en/latest/index.html).
- This can be performed in two steps:
    1. Create a salmon index for the transcriptome
    2. Quantification

## Creating a salmon index <a name="creating_salmon_index"></a>
- To create a *decoy-aware* salmon index requires:
    - genome fasta
    - transcriptome fasta
    - annotation file (GTF)
- A decoy-aware transcriptome is recomended for selective alignment in Salmon from release 0.13.0 onwards.  
- *Why do we need a decoy-aware transcriptome?* When a sequenced fragment originates from an unannotated genomic locus bearing sequence-similarity to an annotated transcript, it can be falsely mapped to the annotated transcript since the relevant genomic sequence is not available to the method. Using MashMap, such sequence-similar decoy regions can be identified and extracted from the genome. The normal Salmon index is then augmented with these decoy sequences, which are handled in a special manner during mapping and alignment scoring, leading to a reduction in false mappings of such a nature.
- Note: check directory, `/tools/salmon/salmonReferences/`, to see whether a salmon index has already been created.
- If not, example code is provided in the next section.


### Example of creating a salmon index <a name="example_salmon_index"></a>
1. Transcriptome fasta was created by concatenating the cDNA and ncRNA fastas from ensembl v97. 
```{bash, eval = F, echo = T}
cd /data/references/fasta/transcriptome/

wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

# Renamed to include version number
gunzip Homo_sapiens.GRCh38.97.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.97.ncrna.fa.gz

# Concatenated to one file
cat Homo_sapiens.GRCh38.97.cdna.all.fa Homo_sapiens.GRCh38.97.ncrna.fa > Homo_sapiens.GRCh38.97.cdna.all.ncrna.fa

# Used wc -l on all three files to check number of lines in concatenate file = number of lines in cdna + ncrna
```
2. Creating a decoy transcriptome file. 
```{bash, eval = F, echo = T}
cd /tools/salmon/salmonReferences/ensembl_v97

bash /tools/salmon/SalmonTools/scripts/generateDecoyTranscriptome.sh \
-j 10 \ 
-b /tools/bedtools2/bin/bedtools \
-m /tools/MashMap/MashMap-2.0/mashmap \
-a /data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf \
-g /data/references/fasta/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa \
-t /data/references/fasta/transcriptome/Homo_sapiens.GRCh38.97.cdna.all.ncrna.fa \
-o /tools/salmon/salmonReferences/ensembl_v97/ 

```
3. Generating the decoy-aware salmon index.
```{bash, eval = F, echo = T}
/tools/salmon/bin/salmon index \
--transcripts /tools/salmon/salmonReferences/ensembl_v97/gentrome.fa \
--index /tools/salmon/salmonReferences/ensembl_v97/Homo_sapiens.GRCh38.97.cdna.all.ncrna_index \
--decoys /tools/salmon/salmonReferences/ensembl_v97/decoy.txt \
--kmerLen 31 # Recommended for > 75 bp reads

```

## Quantification <a name="quantification"></a>

### Determining the library type <a name="library_type"></a>
- To run, Salmon requires that the user specifies the library type. Users then have one of two options:
    1. Allow Salmon to automatically infer the library type by looking at the first few thousand reads of each sample. This is currently the default in the provided script [quantification_Salmon.R](quantification_Salmon.R).
    2. Provide the library type. This can be done in the [quantification_Salmon.R](quantification_Salmon.R) script by calling the `-l/--library_type` flag. A user may know their library type based on the library construction or can use the the `infer_experiment.py` script from `RSeQC` to determine strandedness. For an example command, see below.
    
```{bash, eval = F, echo = T}
/tools/RSeQC-2.6.4/scripts/infer_experiment.py \
-i PDC87_A1B4_GM-T_Aligned.sortedBysamtools.out.bam \
-r /data/references/ensembl/bed/v97/ensembl_GRCh38_v97.bed \
-s 10000000

# Output was majority fraction (> 80%) of reads: 1++,1--,2+-,2-+.
# According to Salmon, this corresponds to ISF.

```
- If in doubt how to interpret the output of `infer_experiment.py`, see [here](https://training.galaxyproject.org/archive/2019-05-01/topics/transcriptomics/tutorials/ref-based/tutorial.html).


### Running Salmon quantification <a name="running_salmon_quant"></a>
Example command shown below.

```{bash, eval = F, echo = T}
nohup Rscript \
/home/rreynolds/packages/RNAseqProcessing/quantification/quantification_Salmon.R \
/data/RNAseq_PD/tissue_polyA_samples/QC/fastp \
/tools/salmon/salmonReferences/ensembl_v97/Homo_sapiens.GRCh38.97.cdna.all.ncrna_index \
/data/RNAseq_PD/tissue_polyA_samples/salmon_quant \
--sample_prefix=NM...._ \
--sample_suffix=_S.* \
--library_type=ISF \
&>/home/rreynolds/projects/Aim2_PDsequencing_wd/Aim2_PDsequencing/nohup_logs/PD_tissue_polyA_salmon_quant.log&
```

### Loading the data post-quantification <a name="loading_post_quant_data"></a>
- Data can be loaded into R using the packages `tximport` and `DESeq2`. As we have R version 3.4 (and thus an older version of Bioconductor), these must be installed using the following commands:
```{R, eval = F, echo = T}
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("tximport")
```
- To import data you need:
    1. A dataframe detailing all of your sample details e.g. age, sex, RIN, etc.
    2. A named vector, wherein names refer to sample IDs, and values of the vector refer to file paths where Salmon `quant.sf` files are stored. 
    3. A two-column dataframe linking transcript ids to gene ids. This is required for methods (such as Salmon) that provide transcript-level estimates only. **N.B. The column names are not relevant, but the order of transcript id and gene id must be used.** 

#### Creating a named vector <a name="creating_named_vector"></a>
If [quantification_Salmon.R](quantification_Salmon.R) used, then `quant.sf` files will be stored in a main `salmon_quant/` directory, which contains sub-directories, each of which carries a sample name. Thus, the following example code can be adapted to create a named vector:
```{R, eval = F, echo = T}
# Create dataframe of file paths and sample names
file_df <- tibble(file_paths = list.files(path = "/path_to_files/salmon_quant", 
                                          recursive = T, 
                                          pattern = "quant.sf",full.names = T),
                  sample_name = list.files(path = "/path_to_files/salmon_quant", 
                                           recursive = T, 
                                           pattern = "quant.sf",full.names = T) %>% 
                    str_replace("/quant.sf", "") %>% 
                    str_replace("/.*/", "")) 

# Create named vector                    
files <- file_df$file_paths
names(files) <- file_df$sample_name
```

#### Creating a transcript-to-gene map <a name="transcript_to_gene_map"></a>
- The easiest way to generate this is to use a `txdb` generated from the reference transcriptome used in your Salmon quantification. 
- If ensembl was used, please check the directory `/data/references/ensembl/txdb_sqlite/` to see if this is already available before creating one.
- If unavailable, a txdb object can be created in a number of ways.
    - It can be prepared from a `.GTF` file, using the `makeTxDbFromGFF` function in the `GenomicFeatures` package.
    - If it is an ensembl transcriptome, it can be generated using the following code:
```{R, eval = F, echo = T}
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                    release=97,
                    server="ensembldb.ensembl.org")
saveDb(txdb,
       file = "/data/references/ensembl/txdb_sqlite/v97/Homo_sapiens.GRCh38.97.sqlite")
```
- Once a txdb object is available, the following code can be used to import the data.
```{R, eval = F, echo = T}
tx2gene <- ensembldb::select(txdb, 
                             keys = keys(txdb, keytype = "TXNAME"), 
                             columns = c("GENEID", "TXNAME"),
                             keytype = "TXNAME")

# Loading in salmon data
txi <- tximport::tximport(files = files,
                          type = c("salmon"),
                          tx2gene = tx2gene, 
                          ignoreTxVersion = TRUE, # splits tx id on the "." character to remove version info for easier matching with tx id in tx2gene
                          dropInfReps = TRUE)
```

### Post-quantification QC checks <a name="QC_checks"></a>
- Once quantification has been performed and data loaded into R, the user can now perform basic quality-control checks, including:
    - Sex checks -- using gene expression for *XIST* (female-specific, expressed for X-chromosome inactivation) and *DDX3Y* (located on Y-chromosome) determine whether gene expression patterns of sex-specific genes match the sample's assigned sex.
    - Check library sizes and count distributions by plotting. Consider colour by various factors (e.g. RIN) to see whether this explains distributions.
    - Principal component (PCA) plot for outlier analysis. If experiment is well controlled and has worked as expected, the greatest sources of variation should arise from the treatment/groups of interest. 
- For a tutorial on the latter two, please refer to this [link](http://sbc.shef.ac.uk/workshops/2018-07-10-rna-seq/rna-seq-preprocessing.nb.html#quality_control) for more details. Code provided can be adapted to fit with your experimental set-up.
