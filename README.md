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
	samtools sort -@ 8 -o 4_H460_2G_1.bam 4_H460_4G_1.sam
	samtools sort -@ 8 -o 5_H460_2G_2.bam 5_H460_4G_2.sam
	samtools sort -@ 8 -o 6_H460_2G_3.bam 6_H460_4G_3.sam
	samtools sort -@ 8 -o 7_hiNeuroS-TRAIL_1.bam 7_hiNeuroS-TRAIL_1.sam
	samtools sort -@ 8 -o 8_hiNeuroS-TRAIL_2.bam 8_hiNeuroS-TRAIL_2.sam
	samtools sort -@ 8 -o 9_hiNeuroS-TRAIL_3.bam 9_hiNeuroS-TRAIL_3.sam
	samtools sort -@ 8 -o 10_hiNeuroS-TRAIL_2G_1.bam 10_hiNeuroS-TRAIL_2G_1.sam
	samtools sort -@ 8 -o 11_hiNeuroS-TRAIL_2G_2.bam 11_hiNeuroS-TRAIL_2G_2.sam
	samtools sort -@ 8 -o 12_hiNeuroS-TRAIL_2G_3.bam 12_hiNeuroS-TRAIL_2G_3.sam
	
	


## index the sam files to generate .bai files

make sure that the only files in the directory are the sam and bam/bai files, then:
	
	find *.bam -exec echo samtools index {} \; ] sh

which essentially converts all sam files automatically, accomplishes the command as:

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


## Use Stringtie to generate expression estimates (FPKM and TPM) from the SAM/BAM files generated by HISAT2
this takes the bam file, add expression estimates FPKM and TPM using the reference annotation GTF file as a guide, and gives an annotated gtf file for the sample.
since the forward read of the resulting sequencing data represents the anti-sense strand, the "-RF" orientation flag should be used. 

	mkdir -p expression/stringtie/ref_only/
	
	cd expression/stringtie/ref_only/
	
	stringtie --rf -p 8 -G $hingtgen/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 1_H460_1/transcripts.gtf -A 1_H460_1/gene_abundances.tsv $hingtgen/alignments/1_H460_1.bam

’–rf’ tells StringTie that our data is stranded and to use the correct strand specific mode
‘-p 8’ tells StringTie to use eight CPUs
‘-G ' reference annotation to use for guiding the assembly process (GTF/GFF3)
‘-e’ only estimate the abundance of given reference transcripts (requires -G)
‘-B’ enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
‘-o’ output path/file name for the assembled transcripts GTF (default: stdout)
‘-A’ output path/file name for gene abundance estimates


View transcript records only and improve formatting

	grep -v "^#" 1_H460_1/transcripts.gtf | grep -w "transcript" | column -t | less -S
	

Limit the view to transcript records and their expression values (FPKM and TPM values)

	awk '{if ($3=="transcript") print}' 1_H460_1/transcripts.gtf | cut -f 1,4,9 | less -S


