# This is the niccolo sample overall pipeline - the framework that will probably be similar for all samples 

# pre alignment QC and trimming reads for - works directly on the fastq files, will output a multiqc 
nohup Rscript \
/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/QC_RNAseq_samples/pre_alignment_QC_trimmomatic_fastQC_multiQC.R \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/FASTQ/C8VP1ANXX/ \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC \
nicc_samples \
&>/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/QC_RNAseq_samples/nohup_logs/nicc_samples.log&

# align the trimmed FASTQ files using STAR 2.pass 
nohup Rscript \
/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/alignment/STAR_alignment.R \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC/trimmomatic/ \
/data/STAR_data/genome_index_hg38 \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR \
&>/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/alignment/nohup_logs/nicc_samples_STAR.log&

# sort and index the bam files using samtools required for some downstream QC
# STAR seems to sort these however, this sorting does not allow for some QC modules such as geneCoverage.py from RSeQC
samtools sort -m 1000000000 NIAA_Aligned.sortedByCoord.out.bam -o NIAA_Aligned.sortedBysamtools.out.bam
samtools sort -m 1000000000 NISA_Aligned.sortedByCoord.out.bam -o NISA_Aligned.sortedBysamtools.out.bam
samtools sort -m 1000000000 NICP_Aligned.sortedByCoord.out.bam -o NICP_Aligned.sortedBysamtools.out.bam
samtools index NIAA_Aligned.sortedBysamtools.out.bam
samtools index NISA_Aligned.sortedBysamtools.out.bam
samtools index NICP_Aligned.sortedBysamtools.out.bam

# remove the trimmomatic fastqs and the STAR aligned .bam files to avoid taking up redundant stace on server
rm -r /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC/trimmomatic
rm /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/*sortedByCoord.out.bam

# run post alignment QC on the bam files
nohup Rscript \
/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/QC_RNAseq_samples/post_alignment_QC_RSeQC.R \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/ \
/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC/ \
/data/references/ensembl/bed/v95/ensembl_GRCh38_v95.bed \
100 \
&>/home/dzhang/projects/transcriptomics_diagnostics_wd/transcriptomics_diagnostics/QC_RNAseq_samples/nohup_logs/nicc_samples_post_alignment_QC.log&

multiqc /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/ -o /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/QC/ \
-n multiqc_niccolo_all_QC

# check the experiment type 
/home/dzhang/.local/bin/infer_experiment.py -i /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/NIAA_Aligned.sortedBysamtools.out.bam \
-r /data/references/ensembl/bed/v95/ensembl_GRCh38_v95.bed
/home/dzhang/.local/bin/infer_experiment.py -i /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/NICP_Aligned.sortedBysamtools.out.bam \
-r /data/references/ensembl/bed/v95/ensembl_GRCh38_v95.bed

# convert the bam's to bw's for derfinder
# ran this is parrellel manually using screen 
/home/dzhang/.local/bin/bam2wig.py -i /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/NIAA_Aligned.sortedBysamtools.out.bam \
-s /data/references/chr_info/hg38.chrom.sizes.ens \
-o /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/bw/NIAA_Aligned.sortedBysamtools -u \
-d '1+-,1-+,2++,2--'
/home/dzhang/.local/bin/bam2wig.py -i /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/NICP_Aligned.sortedBysamtools.out.bam \
-s /data/references/chr_info/hg38.chrom.sizes.ens \
-o /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/bw/NICP_Aligned.sortedBysamtools -u \
-d '1+-,1-+,2++,2--'
/home/dzhang/.local/bin/bam2wig.py -i /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/NISA_Aligned.sortedBysamtools.out.bam \
-s /data/references/chr_info/hg38.chrom.sizes.ens \
-o /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/bw/NISA_Aligned.sortedBysamtools -u \
-d '1+-,1-+,2++,2--'

# leafcutter 
# converting BAMs to .junc files
for bamfile in `ls /data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/STAR/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh /tools/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc
done


