## Comparison of trimmomatic and fastp
Files used for this project were provided by David Zhang:
1. ```/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/FASTQ/C8VP1ANXX/NIAA_*.fastq.gz```
2. ```/data/RNAseq_diagnostics/niccolo_X_linked_dystonia/180420-140039/FASTQ/C8VP1ANXX/NICP_*.fastq.gz```
3. ```/data/RNAseq_diagnostics/mito_samples/postive_controls/2018-09-11/NG-17299_S1483_*.fastq.gz```

Scripts were also adapted from David, and have been written such that they can be generalised to any samples. Depending on whether you want to use fastp or trimmomatic please use:
- Trimmomatic: [prealignmentQC_trimmomatic_fastQC_multiQC.R](prealignmentQC_trimmomatic_fastQC_multiQC.R)
- fastp: [prealignmentQC_fastp_noperreadcutting.R](prealignmentQC_fastp_noperreadcutting.R)
- fastp with trimmomatic settings: [prealignmentQC_fastp.R](prealignmentQC_fastp.R)

A comparison and explanation of these various tools, in addition to a general description of the workflow of QC can be seen using this [link](prealignmentQC_comparison_workflow.html).
