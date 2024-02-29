# Sequencing Data Analysis Toolkit
1. [How to calculate number of reads?](#question1)
2. [How to calculate the average read depth after aligning reads to a reference genome?](#question2)
3. [How to fetch mapped and unmapped reads from BAM file using samtools?](#question3)



This GitHub repository provides a comprehensive collection of commands and scripts commonly used while analyzing sequencing data.

## How to calculate number of reads? <a name="question1"></a>

In a Fastq file, each read is represented by four lines, with each read beginning with the "@" symbol. This command functions by counting the occurrences of "@" at the start of each read.

**For a single file:**
```bash
grep -c "^@" input.fastq
```
**For multiple files:**
```bash
grep -c "^@" *.fastq > read_count.txt
```
**For zipped files:**
```bash
zcat *fastq.gz | grep -c "^@" > read_count.txt
```
## How to calculate the average read depth after aligning reads to a reference genome? <a name="question2"></a>
Depth generally refers to the number of times a particular nucleotide in the genome is read during the sequencing process. It reflects the number of reads aligned to a specific position in the genome. Evaluating read depth helps determine if there are a sufficient number of reads aligned over a given base position, which is crucial for our analysis. Whereas coverage means to check if the reads are aligned over the genome by checking the percentage.
This is a common way to evaluate the quality of samples by checking read depth after alignment.

**Step 1: Convert SAM file to BAM file**
```bash
samtools view -bS --threads <number_of_threads> input.sam > output.bam
```
**Step 2: Sort BAM file based on position** 
```bash
samtools sort input.bam -o output_sorted.bam -@ <number_of_threads>
```
**Step 3: Calculating Average Depth**  
There are two tools either “samtools depth” or “Mosdepth”
```bash
samtools depth -a output_sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

OR 

Index Sorted BAM file and use mosdepth to calculate average read depth
```bash
samtools index input_sorted.bam -@ <number_threads>
mosdepth  sample_depth input_sorted.bam
awk '$1=="total"{print $4;}'  sample_depth.mosdepth.summary.txt
```
## How to fetch mapped and unmapped reads from BAM file using samtools? <a name="question3"></a>
To fetch mapped reads:
```bash
samtools view -b -F 4 input.bam > output.bam
```
“-F” flag to fetch mapped reads & “-b” to produce output BAM file
To fetch unmapped reads:
```bash
samtools view -b -f 4 input.bam > output.bam
```
“-f” flag to fetch unmapped reads

