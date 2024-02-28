# Sequencing Data Analysis Toolkit

This GitHub repository provides a comprehensive collection of commands and scripts commonly used while analyzing sequencing data.

## How to Calculate Number of Reads?

In a Fastq file, each read is represented by four lines, with each read beginning with the "@" symbol. This command functions by counting the occurrences of "@" at the start of each read.

**For a single file:**
```bash
grep -c "^@" <input_fastq_file>
```
**For multiple files:**
```bash
grep -c "^@" *.fastq > read_count.txt
```
**For zipped files:**
```bash
zcat *fastq.gz | grep -c "^@" > read_count.txt
```
## How to calculate the average read depth after aligning reads to a reference genome?
This is a common way to evaluate quality of samples by checking read depth after alignment
Depth generally refers to the number of times a particular nucleotide in the genome is read during the sequencing process. It reflects the number of reads aligned to a specific position in the genome. Evaluating read depth helps determine if there are a sufficient number of reads aligned over a given base position, which is crucial for our analysis. Whereas coverage means to check if the reads are aligned over the genome by checking the percentage.
**Step 1: Convert SAM file to BAM file**
```bash
samtools view -bS --threads <number_of_threads> <input_SAM_file> > <output_bam_file>
```
**Step 2: Sort BAM file based on position** 
```bash
samtools sort <input_BAM> -o <output_sorted_BAM> -@ <number_of_threads>
```
**Step 3: Calculating Average Depth**  
There are two tools either “samtools depth” or “Mosdepth”
```bash
samtools depth -a <output_sorted_bam> |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
OR 
Index Sorted BAM file and use mosdepth to calculate average read depth
```bash
samtools index <input_sorted_bam> -@ <number_threads>
mosdepth  sample_depth <input_sorted_bam>
awk '$1=="total"{print $4;}'  sample_depth.mosdepth.summary.txt
```
