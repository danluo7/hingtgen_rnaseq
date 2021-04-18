# hingtgen_rnaseq
bulk mRNA sequencing 

2021-04-15
downloaded 12 fastq files from Alison Mercer-Smith.

    fastqc *.fastq.gz    / * = wildcard. checked quality of sequencing reads

    zcat MS001_1_S1_R1_001.fastq.gz | grep -P "^\@NS  | wc -l  / checking number of reads by searching for read name prefix and pipe into wc to do line count (-l). Got about 28 million reads per sample.
    
    zcat MS001_1_S1_R1_001.fastq.gz | head -n 8   /take note of the barcode used, will need this to generate the hisat2 bam files later.

Sources for obtaining gene annotation files formatted for HISAT2: HISAT2 Precomputed Genome Index (used by Andrew in Hingtgen's lab) available from their FTP site ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/

    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz
    tar -xzvf hg38.tar.gz

These 8 files together constitute the index: they are all that is needed to align reads to that reference. The original sequence FASTA files are no longer used by HISAT2 once the index is built. So the .fa files and .ht2 files serve the same purpose for HISAT2. If you have a .fa file, then feed this genome.fa fasta file to hisat2-build, and index_name will be the base name of the index files (*.ht2).* Tried to download a genome.fa file to build it, hisat2 killed it after it checked amount of memory available. Why does HISAT2 need 160-200gb to index the human genome? "If you use --snp, --ss, and/or --exon, hisat2-build will need about 200GB RAM for the human genome size as index building involves a graph construction. Otherwise, you will be able to build an index on your desktop with 8GB RAM."

Next is to also download the reference genome from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/. File used by Andrew is: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz 
    gzip -d hg38.ncbiRefSeq.gtf.gz

next is to trimming the FASTQ files if the adaptor sequences have not been trimmed off. check this in fastqc (last output graph shows adaptor sequence content). if need trimming:


    mkdir FASTQ_trimmed
    wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa  /this is to get the adaptor sequences from illumina

Use flexbar to remove illumina adapter sequences.The left side of reads is kept if long enough. The overlap of an adapter and read must have at least length 3 with at most 10% errors in default settings. since this is single-end:

    flexbar -r $hingtgen/FASTQs/MS001_1_S1_R1_001.fastq.gz -a $hingtgen/FASTQs_trimmed/illumina_multiplex.fa -ao 3 -ae 0.1 --target $hingtgen/FASTQs_trimmed/H460_1

next is alignment using hisat2 (which suceeded hisat and tophat)

setting the Read Group info:

a) Sample and library tags. Can be autmatically detected from current sample naming scheme:<Sample.ID><Index.Sequence><Lane.ID><Set.number>.fastq

which for us is: MS001_1_S1_R1_001.fastq.gz
sample.ID: MS001
index.sequence: 1
lane.id: S1
set.number: R1


SM = <Sample.ID>
LB = <Sample.ID>_<Index.Sequence>


b) ID and PU (to enable merging replictes)

ID = Read group identifier = {FLOWCELL_BARCODE}.{LANE}
PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.

Also can be identified from the name of a sequence read in the Fastq file:
@(instrument id)
:(run number)
:(flowcell ID)
:(lane)
:(tile)
:(x_pos)
:(y_pos) (read)
:(is filtered)
:(control number)
:(index sequence)FLOWCELL_BARCODE = @(instrument id)
:(run number)
:(flowcell ID)


    hisat2 -p 8 --rg-id=H460_1 --rg SM:MS001 --rg LB:MS001_1 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.AACCAGAG -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_1_S1_R1_001.fastq.gz -S $hingtgen/alignments/1_H460_1.sam

a few things to note: 1. the reference genome is the genome.1-6.ht2 files, so when pointing to them, use directory/genome without the file extensions. 2) since this is single end read, no need to use -1 option in front of the input file.


    hisat2 -p 8 --rg-id=H460_2 --rg SM:MS001 --rg LB:MS001_2 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.GTCAGTCA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_2_S2_R1_001.fastq.gz -S $hingtgen/alignments/2_H460_2.sam
    
    hisat2 -p 8 --rg-id=H460_3 --rg SM:MS001 --rg LB:MS001_3 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CCTTCCAT -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_2_S2_R1_001.fastq -S $hingtgen/alignments/3_H460_3.sam

