# hingtgen_rnaseq project: bulk mRNA sequencing 

2021-04-15
downloaded 12 fastq files from Alison Mercer-Smith to Floyd lab server. Path: All_Staff/Raw Data/ARM-mRNAseq_Hingtgen

## Check quality of sequencing reads
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





### next is to trim the FASTQ files if the adaptor sequences have not been trimmed off. 
check this in fastqc (last output graph shows adaptor sequence content). if need trimming:


    mkdir FASTQ_trimmed
    wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa  /this is to get the adaptor sequences from illumina

Use flexbar to remove illumina adapter sequences.The left side of reads is kept if long enough. The overlap of an adapter and read must have at least length 3 with at most 10% errors in default settings. since this is single-end:

    flexbar -r $hingtgen/FASTQs/MS001_1_S1_R1_001.fastq.gz -a $hingtgen/FASTQs_trimmed/illumina_multiplex.fa -ao 3 -ae 0.1 --target $hingtgen/FASTQs_trimmed/H460_1



### decoding samples with key from Alison:
Here’s the key with the RNA concentrations I submitted in ng/uL:
				1 - H460, 140
				2 - H460, 99
				3 - H460, 98
				4 - 2 Gy H460, 100
				5 - 2 Gy H460, 98
				6 - 2 Gy H460, 94
				7 - hiNeuroS-TRAIL, 95
				8 - hiNeuroS-TRAIL, 84
				9 - hiNeuroS-TRAIL, 84
				10 - 2 Gy hiNeuroS-TRAIL, 86
				11 - 2 Gy hiNeuroS-TRAIL, 93
				12 - 2 Gy hiNeuroS-TRAIL, 89


### next is alignment using hisat2 (which suceeded hisat and tophat)

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


a few things to note: 1) the reference genome is the genome.1-6.ht2 files, so when pointing to them, use directory/genome without the file extensions. 2) since this is single end read, no need to use -1 option in front of the input file.

    hisat2 -p 8 --rg-id=H460_1 --rg SM:MS001 --rg LB:MS001_1 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.AACCAGAG -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_1_S1_R1_001.fastq.gz -S $hingtgen/alignments/1_H460_1.sam
    hisat2 -p 8 --rg-id=H460_2 --rg SM:MS001 --rg LB:MS001_2 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.GTCAGTCA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_2_S2_R1_001.fastq -S $hingtgen/alignments/2_H460_2.sam
    hisat2 -p 8 --rg-id=H460_3 --rg SM:MS001 --rg LB:MS001_3 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CCTTCCAT -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_3_S3_R1_001.fastq -S $hingtgen/alignments/3_H460_3.sam

    hisat2 -p 8 --rg-id=H460_4G_1 --rg SM:MS001 --rg LB:MS001_4 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.AGGAACAC -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_4_S4_R1_001.fastq -S $hingtgen/alignments/4_H460_2G_1.sam
    hisat2 -p 8 --rg-id=H460_4G_2 --rg SM:MS001 --rg LB:MS001_5 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CTTACAGC -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_5_S5_R1_001.fastq -S $hingtgen/alignments/5_H460_2G_2.sam
    hisat2 -p 8 --rg-id=H460_4G_3 --rg SM:MS001 --rg LB:MS001_6 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.TACCTGCA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_6_S6_R1_001.fastq -S $hingtgen/alignments/6_H460_2G_3.sam
     
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_1 --rg SM:MS001 --rg LB:MS001_7 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.AGACGCTA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_7_S7_R1_001.fastq -S $hingtgen/alignments/7_hiNeuroS-TRAIL_1.sam  
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_2 --rg SM:MS001 --rg LB:MS001_8 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CAACACAG -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_8_S8_R1_001.fastq -S $hingtgen/alignments/8_hiNeuroS-TRAIL_2.sam   
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_3 --rg SM:MS001 --rg LB:MS001_9 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.GTACCACA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_9_S9_R1_001.fastq -S $hingtgen/alignments/9_hiNeuroS-TRAIL_3.sam
        
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_2G_1 --rg SM:MS001 --rg LB:MS001_10 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CGAATACG -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_10_S10_R1_001.fastq -S $hingtgen/alignments/10_hiNeuroS-TRAIL_2G_1.sam    
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_2G_2 --rg SM:MS001 --rg LB:MS001_11 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.GTCCTTGA -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_11_S11_R1_001.fastq -S $hingtgen/alignments/11_hiNeuroS-TRAIL_2G_2.sam     
    hisat2 -p 8 --rg-id=hiNeuroS-TRAIL_2G_3 --rg SM:MS001 --rg LB:MS001_12 --rg PL:ILLUMINA --rg PU:HW3MNBGXH.1.CAGTGCTT -x $hingtgen/RNA_REF_FA/hg38/genome --dta --rna-strandness RF $hingtgen/FASTQs/MS001_12_S12_R1_001.fastq -S $hingtgen/alignments/12_hiNeuroS-TRAIL_2G_3.sam


## convert SAM files to BAM files (saves tons of space) and sort by aligned position

	cd $hingtgen/alignments
	
	samtools sort -@ 8 -o 1_H460_1.bam 1_H460_1.sam
	samtools sort -@ 8 -o 2_H460_2.bam 2_H460_2.sam
	samtools sort -@ 8 -o 3_H460_3.bam 3_H460_3.sam
	
	samtools sort -@ 8 -o 4_H460_2G_1.bam 4_H460_2G_1.sam
	samtools sort -@ 8 -o 5_H460_2G_2.bam 5_H460_2G_2.sam
	samtools sort -@ 8 -o 6_H460_2G_3.bam 6_H460_2G_3.sam
	
	samtools sort -@ 8 -o 7_hiNeuroS-TRAIL_1.bam 7_hiNeuroS-TRAIL_1.sam
	samtools sort -@ 8 -o 8_hiNeuroS-TRAIL_2.bam 8_hiNeuroS-TRAIL_2.sam
	samtools sort -@ 8 -o 9_hiNeuroS-TRAIL_3.bam 9_hiNeuroS-TRAIL_3.sam
	
	samtools sort -@ 8 -o 10_hiNeuroS-TRAIL_2G_1.bam 10_hiNeuroS-TRAIL_2G_1.sam
	samtools sort -@ 8 -o 11_hiNeuroS-TRAIL_2G_2.bam 11_hiNeuroS-TRAIL_2G_2.sam
	samtools sort -@ 8 -o 12_hiNeuroS-TRAIL_2G_3.bam 12_hiNeuroS-TRAIL_2G_3.sam
	
	


## index the sam files to generate .bai files

make sure that the only files in the directory are the sam and bam/bai files, then:
	
	find *.bam -exec echo samtools index {} \; | sh

which essentially converts all sam files automatically, accomplishes the same commands as:

	samtools index 1_H460_1.bam
	samtools index 2_H460_2.bam...


    
## visualize with IGV
start IGV from shell script and navigate to folder with the .bam and bai files
Notice that all of the reads are pointing from 5'-3', another indication of this being a SE (single-end) read seq, with all reads aligning to the sense strand, and therefore the resulting data represents the "anti-sense strand". 

## alignment QC
taking what was aligned by HISAT2 and QC them

	samtools view -H 1_H460_1.bam
	samtools view 1_H460_1.bam | head
	samtools view 1_H460_1.bam | head | column -t | less -S    / this puts it in nicer looking columns 
exit with "q"
	
then filter OUT all reads that are unmapped, mate is unmapped, and not primary alignment. flags that we want. flags numbers are in this webpage: http://broadinstitute.github.io/picard/explain-flags.html

	samtools view -F 260 1_H460_1.bam | head | column -t | less -S

## use samtoools flagstat to geta basic sumary of an alignment 
for example percent of unmapped reads

	mkdir flagstat
	samtools flagstat 1_H460_1.bam > flagstat/1_H460_1.bam.flagstat

view the resulting flagstat: 
	
	cat flagstat/1_H460_1.bam.flagstat
	
output: 
45179552 + 0 in total (QC-passed reads + QC-failed reads)
6992528 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
42501017 + 0 mapped (94.07% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)


## Use Stringtie to generate expression estimates (FPKM and TPM) from the SAM/BAM files generated by HISAT2 and samtools
this takes the bam file, add expression estimates FPKM and TPM using the reference annotation GTF file as a guide, and gives an annotated gtf file for the sample.
since the forward read of the resulting sequencing data represents the anti-sense strand, the "-RF" orientation flag should be used. 

	mkdir -p expression/stringtie/ref_only/
	
	cd expression/stringtie/ref_only/
	
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 1_H460_1/transcripts.gtf -A 1_H460_1/gene_abundances.tsv $hingtgen/alignments/1_H460_1.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 2_H460_2/transcripts.gtf -A 2_H460_2/gene_abundances.tsv $hingtgen/alignments/2_H460_2.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 3_H460_3/transcripts.gtf -A 3_H460_3/gene_abundances.tsv $hingtgen/alignments/3_H460_3.bam
	
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 4_H460_2G_1/transcripts.gtf -A 4_H460_2G_1/gene_abundances.tsv $hingtgen/alignments/4_H460_2G_1.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 5_H460_2G_2/transcripts.gtf -A 5_H460_2G_2/gene_abundances.tsv $hingtgen/alignments/5_H460_2G_2.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 6_H460_2G_3/transcripts.gtf -A 6_H460_2G_3/gene_abundances.tsv $hingtgen/alignments/6_H460_2G_3.bam
	
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 7_hiNeuroS-TRAIL_1/transcripts.gtf -A 7_hiNeuroS-TRAIL_1/gene_abundances.tsv $hingtgen/alignments/7_hiNeuroS-TRAIL_1.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 8_hiNeuroS-TRAIL_2/transcripts.gtf -A 8_hiNeuroS-TRAIL_2/gene_abundances.tsv $hingtgen/alignments/8_hiNeuroS-TRAIL_2.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 9_hiNeuroS-TRAIL_3/transcripts.gtf -A 9_hiNeuroS-TRAIL_3/gene_abundances.tsv $hingtgen/alignments/9_hiNeuroS-TRAIL_3.bam
	
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 10_hiNeuroS-TRAIL_2G_1/transcripts.gtf -A 10_hiNeuroS-TRAIL_2G_1/gene_abundances.tsv $hingtgen/alignments/10_hiNeuroS-TRAIL_2G_1.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 11_hiNeuroS-TRAIL_2G_2/transcripts.gtf -A 11_hiNeuroS-TRAIL_2G_2/gene_abundances.tsv $hingtgen/alignments/11_hiNeuroS-TRAIL_2G_2.bam
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 12_hiNeuroS-TRAIL_2G_3/transcripts.gtf -A 12_hiNeuroS-TRAIL_2G_3/gene_abundances.tsv $hingtgen/alignments/12_hiNeuroS-TRAIL_2G_3.bam
	
	

’–rf’ tells StringTie that our data is stranded and to use the correct strand specific mode
‘-p 8’ tells StringTie to use eight CPUs
‘-G ' reference annotation to use for guiding the assembly process (GTF/GFF3)
‘-e’ only estimate the abundance of given reference transcripts (requires -G)
‘-B’ enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
‘-o’ output path/file name for the assembled transcripts GTF (default: stdout)
‘-A’ output path/file name for gene abundance estimates


view raw output from stringtie

	less -S 1_H460_1/transcripts.gtf
	
	
View transcript records only (column -t makes formatting easier to view). scroll right to see the ### 3 major expression matrices "cov, fpkm, and tpm". ### 
Cov, or coverage, is average per base covereage over the genomic segment. basically, the coverage of each transcript. it's  un-normalized unlike fpkm or tpm. 

	grep -v "^#" 1_H460_1/transcripts.gtf | grep -w "transcript" | column -t | less -S

	

Limit the view only to transcript records and their expression estimates (FPKM and TPM values). the if statement {if ($3=="transcript") print} means only print the line if the third column is a transcript. then only return field 1, 4, 9 which are the fpkm nad tpm values.

	awk '{if ($3=="transcript") print}' 1_H460_1/transcripts.gtf | cut -f 1,4,9 | less -S

Could also view gene and trnscript level expression values in the two files generated by stringtie, a transcript.gtf file (transcript level abundance) and a gene_abundances.tsv (a gene level abundance):

	column -t 1_H460_1/t_data.ctab | less -S
	less -S -x20 1_H460_1/gene_abundances.tsv



(still figuring out the perl script)





## Parallel to Stringtie, can also run htseq-count on alignments to produce raw counts instead of FPKM/TPM values (what Stringtie outputs) for differential expression analysis ##

htseq-count basic usage:

	htseq-count [options] <sam_file> <gff_file>

Extra options specified below:

    ’–format’ specify the input file format one of BAM or SAM. Since we have BAM format files, select ‘bam’ for this option.
    ’–order’ provide the expected sort order of the input file. Previously we generated position sorted BAM files so use ‘pos’.
    ’–mode’ determines how to deal with reads that overlap more than one feature. We believe the ‘intersection-strict’ mode is best.
    ’–stranded’ specifies whether data is stranded or not. The TruSeq strand-specific RNA libraries suggest the ‘reverse’ option for this parameter.
    ’–minaqual’ will skip all reads with alignment quality lower than the given minimum value
    ’–type’ specifies the feature type (3rd column in GFF file) to be used. (default, suitable for RNA-Seq and Ensembl GTF files: exon)
    ’–idattr’ The feature ID used to identify the counts in the output table. The default, suitable for RNA-SEq and Ensembl GTF files, is gene_id.
	

	mkdir -p expression/htseq_counts
	cd expression/htseq_counts

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/1_H460_1.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 1_H460_1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/2_H460_2.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 2_H460_2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/3_H460_3.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 3_H460_3.tsv
	
	
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/4_H460_2G_1.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 4_H460_2G_1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/5_H460_2G_2.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 5_H460_2G_2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/6_H460_2G_3.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 6_H460_2G_3.tsv
	
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/7_hiNeuroS-TRAIL_1.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 7_hiNeuroS-TRAIL_1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/8_hiNeuroS-TRAIL_2.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 8_hiNeuroS-TRAIL_2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/9_hiNeuroS-TRAIL_3.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 9_hiNeuroS-TRAIL_3.tsv
	
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/10_hiNeuroS-TRAIL_2G_1.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 10_hiNeuroS-TRAIL_2G_1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/11_hiNeuroS-TRAIL_2G_2.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 11_hiNeuroS-TRAIL_2G_2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $hingtgen/alignments/12_hiNeuroS-TRAIL_2G_3.bam $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 12_hiNeuroS-TRAIL_2G_3.tsv
	
	

Join the results for each replicate together and merge results files into a single matrix for use in edgeR

	join 1_H460_1.tsv 2_H460_2.tsv | join - 3_H460_3.tsv | join - 4_H460_2G_1.tsv | join - 5_H460_2G_2.tsv |join - 6_H460_2G_3.tsv | join - 7_hiNeuroS-TRAIL_1.tsv | join - 8_hiNeuroS-TRAIL_2.tsv | join - 9_hiNeuroS-TRAIL_3.tsv | join - 10_hiNeuroS-TRAIL_2G_1.tsv | join - 11_hiNeuroS-TRAIL_2G_2.tsv | join - 12_hiNeuroS-TRAIL_2G_3.tsv > gene_read_counts_table_all.tsv
	

Creat a simple text file with just the header that will be used for the table:

	echo "GeneID 1_H460_1 2_H460_2 3_H460_3 4_H460_2G_1 5_H460_2G_2 6_H460_2G_3 7_hiNeuroS-TRAIL_1 8_hiNeuroS-TRAIL_2 9_hiNeuroS-TRAIL_3 10_hiNeuroS-TRAIL_2G_1 11_hiNeuroS-TRAIL_2G_2 12_hiNeuroS-TRAIL_2G_3" > header.txt
	
Clean up a bit more, add a header, reformat the result as a tab delimited file.
note: grep -v "__" is being used to filter out the summary lines at the end of the files that ht-seq count gives to summarize reads that had no feature, were ambiguous, did not align at all, did not align due to poor alignment quality, or the alignment was not unique.

awk -v OFS="\t" '$1=$1' is using awk to replace the single space characters that were in the concatenated version of our header.txt and gene_read_counts_table_all.tsv with a tab character. -v is used to reset the variable OFS, which stands for Output Field Separator. By default, this is a single space. By specifying OFS="\t", we are telling awk to replace the single space with a tab. The '$1=$1' tells awk to reevaluate the input using the new output variable

	cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > gene_read_counts_table_all_final.tsv

	rm -f gene_read_counts_table_all.tsv header.txt

	head gene_read_counts_table_all_final.tsv | column -t

Notice that only about 20k genes have considerable number of reads mapped to them.




## Use Ballgown in R for differential expression (DE) analysis using output from Stringtie ##
## Perform A vs. B comparison, using all replicates, for known (reference only mode) transcripts

(raw counts using htseq output will come later)

	mkdir -p $hingtgen/de/ballgown/ref_only
	cd $hingtgen/de/ballgown/ref_only/



Use printf to create/print a table with ids, type (each sample is a type), and path to the file, as the header. Then n returns a new line.

Bascially, we need a table that needs to look like this to feed into R:

ids		type		path to file
1_H460_1	H460		$hingtgen/expression/stringtie/ref_only/1_H460_1
2_H460_2	H460		$hingtgen/expression/stringtie/ref_only/2_H460_2
...
...

this is how the script should look like (without the enters inbetween each line):

printf "\"ids\",\"type\",\"path\"\

n\"1_H460_1\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/1_H460_1\"\
n\"2_H460_2\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/2_H460_2\"\
n\"3_H460_3\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/3_H460_3\"\

n\"4_H460_2G_1\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/4_H460_2G_1\"\
n\"5_H460_2G_2\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/5_H460_2G_2\"\
n\"6_H460_2G_3\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/6_H460_2G_3\"\

n\"7_hiNeuroS-TRAIL_1\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/7_hiNeuroS-TRAIL_1\"\
n\"8_hiNeuroS-TRAIL_2\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/8_hiNeuroS-TRAIL_2\"\
n\"9_hiNeuroS-TRAIL_3\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/9_hiNeuroS-TRAIL_3\"\

n\"10_hiNeuroS-TRAIL_2G_1\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/10_hiNeuroS-TRAIL_2G_1\"\
n\"11_hiNeuroS-TRAIL_2G_2\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/11_hiNeuroS-TRAIL_2G_2\"\
n\"12_hiNeuroS-TRAIL_2G_3\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/12_hiNeuroS-TRAIL_2G_3\"\




full script:


	printf "\"ids\",\"type\",\"path\"\n\"1_H460_1\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/1_H460_1\"\n\"2_H460_2\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/2_H460_2\"\n\"3_H460_3\",\"H460\",\"$hingtgen/expression/stringtie/ref_only/3_H460_3\"\n\"4_H460_2G_1\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/4_H460_2G_1\"\n\"5_H460_2G_2\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/5_H460_2G_2\"\n\"6_H460_2G_3\",\"H460_2G\",\"$hingtgen/expression/stringtie/ref_only/6_H460_2G_3\"\n\"7_hiNeuroS-TRAIL_1\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/7_hiNeuroS-TRAIL_1\"\n\"8_hiNeuroS-TRAIL_2\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/8_hiNeuroS-TRAIL_2\"\n\"9_hiNeuroS-TRAIL_3\",\"hiNeuroS-TRAIL\",\"$hingtgen/expression/stringtie/ref_only/9_hiNeuroS-TRAIL_3\"\n\"10_hiNeuroS-TRAIL_2G_1\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/10_hiNeuroS-TRAIL_2G_1\"\n\"11_hiNeuroS-TRAIL_2G_2\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/11_hiNeuroS-TRAIL_2G_2\"\n\"12_hiNeuroS-TRAIL_2G_3\",\"hiNeuroS-TRAIL_2G\",\"$hingtgen/expression/stringtie/ref_only/12_hiNeuroS-TRAIL_2G_3\"\n" > all_4_comparisons.csv

start R and load libraries

	R
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	
Load phenotype data from the file we just saved in the current working directory
	
	pheno_data = read.csv("all_4_comparisons.csv")

Load ballgown data structure and save it to a variable "bg"

	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

Display a description of this object

	bg


Load all attributes including gene name, put in a table by using a transcript expression function from the ballgown library (texpr), feed it the ballgown library, all of it. 

	bg_table = texpr(bg, 'all')

Then pull out just the gene names and gene id's. [, 9:10] means pull out all rows and just column 9 and 10. Keep in mind that the ncbi reference genomes gene id and gene names are the same, so column 9 and 10 will be the same. [, x,x] is essentially how to subset data within R. This table is used in later step for DE after merging the results.

	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


Save the ballgown object to a file, R data file, for later use, so later you odn't have to load each sample from their directories.

	save(bg, file='bg.rda')



## Perform differential expression (DE) analysis with no filtering

creat a results object. Use function stattest, feed the ballgown object, tell it whta we want to look at is transcript data, and the covariate that it's gonna perform the DE on is type, which was created earlier in the table using the printf function. type was H460 vs H460_2G etc. Use getFC to get fold change. Use "meas" to use FPKM as measurement of the foldchange


	results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
	
same thing with gene level data

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	
	
	results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))






