# pre alignment QC and trimming reads using fastQC and trimmomatic - works directly on the fastq files, will output a multiqc
nohup Rscript \
/home/rreynolds/projects/Aim2_TestQCandAlignment/QC_RNAseq_samples/prealignmentQC_trimmomatic_fastQC_multiQC.R \
/home/rreynolds/data/RNASeq_TestRunQC \
/home/rreynolds/data/RNASeq_TestRunQC/QC \
testrun_fastQC_trimmomatic \
NG-17299_ \
&>/home/rreynolds/projects/Aim2_TestQCandAlignment/nohup_logs/testrun_fastQC_trimmomatic.log&

# pre alignment QC and trimming reads using fastp with trimmomatic settings - works directly on the fastq files, will output a multiqc
nohup Rscript \
/home/rreynolds/projects/Aim2_TestQCandAlignment/QC_RNAseq_samples/prealignmentQC_fastp.R \
/home/rreynolds/data/RNASeq_TestRunQC \
/home/rreynolds/data/RNASeq_TestRunQC/QC \
testrun_fastp_trimmomatic_settings \
NG-17299_ \
&>/home/rreynolds/projects/Aim2_TestQCandAlignment/nohup_logs/testrun_fastp_trimmomatic_settings.log&

# pre alignment QC and trimming reads using fastp without trimmomatic settings - works directly on the fastq files, will output a multiqc
nohup Rscript \
/home/rreynolds/projects/Aim2_TestQCandAlignment/QC_RNAseq_samples/prealignmentQC_fastp_noperreadcutting.R \
/home/rreynolds/data/RNASeq_TestRunQC \
/home/rreynolds/data/RNASeq_TestRunQC/QC \
testrun_fastp_noperreadcutting \
NG-17299_ \
&>/home/rreynolds/projects/Aim2_TestQCandAlignment/nohup_logs/testrun_fastp_noperreadcutting.log&
