# Running QUACk 

Here is shown an example of how to run QUACk in order to set an exome analysis pipeline. A set of two paired-end library reads have been sequenced for each of two samples

### Command line

    perl /PATH/TO/quack/quack.pl -sample_file test.sample -parameter_file test.param -reference hg19 -project test -threads 1 -fastq_extension fastq -fastq_suffix _001 -extended_bed _extended

**Example files**

- **test.sample**: example sample file
- **test.param**: example parameter file
- **GENCODE_coding.srt.merged.chr1_extended.bed**: extended bed file

**!!!** FASTA file **MUST** be provided by the user
	
---
**Sample file** and **parameter file** are tab-delimited and comma-delimited files, and contain information about sample sequencing features and analysis parameters, respectively. 
