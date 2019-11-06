# Running QUACk 

Here is shown an example of how to run QUACk in order to set up an exome analysis pipeline. A set of two paired-end library reads have been sequenced for each of two samples:

- tp_143h100K
    - tp_143h100K_FCD2CHUACXX_L8_AAAHAAA_R1_001.fastq  
    - tp_143h100K_FCD2HJGACXX_L6_AAAHAAA_R1_001.fastq
    - tp_143h100K_FCD2CHUACXX_L8_AAAHAAA_R2_001.fastq  
    - tp_143h100K_FCD2HJGACXX_L6_AAAHAAA_R2_001.fastq

- tp_148h100K
    - tp_148h100K_DCW97JN1_L1_CGATGT_R1_001.fastq  
    - tp_148h100K_DCW97JN1_L1_CGATGT_R2_001.fastq  
    - tp_148h100K_HISEQ5_L1_CGATGT_R1_001.fastq  
    - tp_148h100K_HISEQ5_L1_CGATGT_R2_001.fastq

### Command line

    perl /PATH/TO/quack/quack.pl -sample_file test.sample -parameter_file test.param -reference hg19 -project test -threads 1 -fastq_extension fastq -fastq_suffix _001 -extended_bed _extended

**Example files**

- **test.sample**: example sample file
- **test.param**: example parameter file
- **GENCODE_coding.srt.merged.chr1_extended.bed**: extended bed file
- **test**: folder containing sample subfolders with test fastQ files

**!!!** FASTA file **MUST** be provided by the user
