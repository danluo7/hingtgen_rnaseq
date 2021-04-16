# hingtgen_rnaseq
bulk mRNA sequencing 

2021-04-15
downloaded 12 fastq files from Alison Mercer-Smith.

fastqc *.fastq.gz    / * = wildcard. checked quality of sequencing reads

zcat MS001_1_S1_R1_001.fastq.gz | grep -P "^\@NS  | wc -l  / checking number of reads by searching for read name prefix and pipe into wc to do line count (-l). Got about 30 million reads per sample.



