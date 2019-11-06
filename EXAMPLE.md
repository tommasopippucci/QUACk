# Running QUACk 

Here is shown an example of how to run QUACk in order to set an exome analysis pipeline. A set of two paired-end library reads have been sequenced for each of two samples

### Command line

    perl /PATH/TO/quack/quack.pl -sample_file test.sample -parameter_file test.param -reference hg19 -project test -threads 1 -fastq_extension fastq -fastq_suffix _001 -extended_bed _extended

**Mandatory options**

        -sample_file <SAMPLE_FILE.sample> 
        -parameter_file <PARAMETER_FILE.param> 
        -reference <hg19|hg38|..> 
        -project <PROJECT_FOLDER>
	
---
**Sample file** and **parameter file** are tab-delimited and comma-delimited files, and contain information about sample sequencing features and analysis parameters, respectively. 
