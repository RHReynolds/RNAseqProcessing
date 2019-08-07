--WORK IN PROGRESS--

## Download + installation 

- Salmon downloaded from https://github.com/alexdobin/STAR/releases. Current version is 0.14.1.
- Placed and unzipped into "/tools/salmon/salmon-latest_linux_x86_64/". These come pre-compiled and the binary executable files are found in "/bin".
- Optional: add the binary executable to PATH.

```{bash download/install salmon, echo = T, tidy = T, eval = F}
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz
tar -xvzf salmon-0.14.1_linux_x86_64.tar.gz
export PATH=$PATH:/tools/salmon/salmon-latest_linux_x86_64/bin/
```

## Running Salmon
- While Salmon can be run using previously aligned .bam files, it requires that these were aligned directly to the transcriptome rather than to the genome. There are ways of converting genome-aligned .bam files; however, given the speed of Salmon mapping, opted to perform quasi-mapping and quantification (otherwise known as mapping-based mode in Salmon documentation - https://salmon.readthedocs.io/en/latest/index.html).
- This can be performed in two steps:
    1. Create a salmon index for the transcriptome
    2. Quantification

## Creating a salmon index
- To create a *decoy-aware* salmon index requires:
    - genome fasta
    - transcriptome fasta
    - annotation file (GTF)
- A decoy-aware transcriptome is recomended for selective alignment in Salmon from release 0.13.0 onwards.  When a sequenced fragment originates from an unannotated genomic locus bearing sequence-similarity to an annotated transcript, it can be falsely mapped to the annotated transcript since the relevant genomic sequence is not available to the method. Using of MashMap, such sequence-similar decoy regions can be identified and extracted from the genome. The normal Salmon index is then augmented with these decoy sequences, which are handled in a special manner during mapping and alignment scoring, leading to a reduction in false mappings of such a nature.
- Note: check directory, `/tools/salmon/salmonReferences/`, to see whether a salmon index has already been created.


### Example of creating a salmon index 
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

## Quantification

### Determining the library type
- To run, Salmon requires that the user specifies the library type. Users then have one of two options:
    1. Allow Salmon to automatically infer the library type by looking at the first few thousand reads of each sample. This is currently the default in the provided script [quantification_Salmon.R](quantification_Salmon.R).
    2. Provide the library type. This can be done in the [quantification_Salmon.R](quantification_Salmon.R) script by calling the `-l/--library_type` flag. A user may know their library type based on the library construction or can use the the `infer_experiment.py` script from `RSeQC` to determine strandedness. For an example command, see below.
    
```{bash, eval = F, echo = T}
/tools/RSeQC-2.6.4/scripts/infer_experiment.py -i PDC87_A1B4_GM-T_Aligned.sortedBysamtools.out.bam -r /data/references/ensembl/bed/v97/ensembl_GRCh38_v97.bed -s 10000000

# Output was majority fraction (> 80%) of reads: 1++,1--,2+-,2-+.
# According to Salmon, this corresponds to ISF.

```
- If in doubt how to interpret the output of `infer_experiment.py`, see:  https://training.galaxyproject.org/archive/2019-05-01/topics/transcriptomics/tutorials/ref-based/tutorial.html


### Running Salmon quantification
