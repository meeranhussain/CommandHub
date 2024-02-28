# Sequencing Data Analysis Toolkit

This GitHub repository provides a comprehensive collection of commands and scripts commonly used for analyzing sequencing data.

## How to Calculate Number of Reads?

There are two methods for calculating the number of reads:

### Method 1: Using Awk

In a Fastq file, each read is represented by four lines, with each read beginning with the "@" symbol. This command functions by tallying the occurrences of "@" at the start of each read.

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
