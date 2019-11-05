# QUACk

QUACk (QUality Alignment and Calling) Perl script generating bash (.sh) scripts to perform three-steps NGS analysis (1. Read quality; 2. Alignment and alignment statistics; 3. Creation of gVCF)

## Usage

**Command**

    quack.pl 

**Mandatory options**

        -sample_sheet <SAMPLE_FILE.sample> 
	-param <PARAMETER_FILE.param> 
	-reference <hg19|hg38|..> 
	-project <PATH/TO/PROJECT/FOLDER> 

**Optional parameters**

        -aligner <default:bwa> 
	-threads <default:1> 
	-fastq_extension <default:fastq.gz> 
	-fastq_suffix <default:no_string>

---

## INPUT files

### *FastQ* files, having the following mandatory filename form:

<Sample_Name>_L<Lane_Number>_<FlowCell_ID>_<Index_Sequence>_R<1(forward)|2(reverse)><fastQ_suffix>.<fastQ_extension>

### *SAMPLE* file, tab-delimited

|Sample | FlowCell | Lane | Index | Enrichment | Target_Set | Library | Platform | Provider|
|---    |---       |---   |---    |---         |---         |---      |---       |---      |
|SAMPLE1|  FCID    |1     |  ACAG |MedExome    |Gencode     |Library1 | Illumina |  Seq    |
|SAMPLE1|  FCID    |2     |  GGAG |MedExome    |Gencode     |Library1 | Illumina |  Seq    |
|SAMPLE2|  FCID    |3     |  AGAG |MedExome    |Gencode     |Library2 | Illumina |  Seq    |
|SAMPLE3|  FCID    |2     |  TCAG |MedExome    |Gencode     |Library3 | Illumina |  Seq    |

---
- Sample: Sample Name
- FlowCell: FlowCell ID
- Lane: FlowCell Lane number
- Index: Index Sequence
- Enrichment: Target enrichment used (if used)
- Target Set: Target set used in the analysis (can be different from Enrichment)
- Library: Name of library performed
- Platform: Sequencing platform used
- Provider: Name of the sequencing provider

---

### *PARAMETER* file, comma-delimited

---
- **reference sequence fasta/fa files** (can be different reference sequences)
  - hg19,/archive/ngsbo/db/hg19/ucsc.hg19.fasta
  - GRCh38,/archive/ngsbo/db/GRCh38/hs38.fa
  - hs37d5,/archive/ngsbo/db/hs37d5/hs37d5.fa
- **wes target bed files folder** (complete folder path where target bed files are stored)
  - target_dir,/archive/ngsbo/db/regions/
- **alignment algorithm commands** (command-line command for aligner (no options, no arguments))
  - bwa,bwa
- **SAM/BAM manipulation and calling algorithm commands** (command-line command for samtools, picard, gatk, ... (no options, no arguments))
  - samtools,samtools
  - picard,picard
  - gatk,gatk
- **(VCF) files with known sites of true variation** (dbSNP (or other) true variation sites to be used)
  - known_sites,/archive/ngsbo/db/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf
- **working directory** (omplete path of directory where tmp and quality files will be stored)
  - working_dir,/work/common/pipeline_wes/
- **fastQ directory** (complete path of directory where fastQ filescan be retrieved)
  - fastq_dir,/archive/ngsbo/fastq/
- **BAM directory** (directory where to store final BAM files)
  - bam_dir,/archive/ngsbo/bam/
- **GVCF directory** (directory where to store gVCF files)
  - gvcf_dir,/archive/ngsbo/gvcf/


---


## Directory architecture

QUACk checks for existence of the working, bam and fastq directories, and generates a directory architecture to store output and temporary files. If the output folder does exist, QUACk generates a folder architecture as follows:

---

    project

        project_analysis [log file]
	
	    bash

                calling [calling.sh]

                alignment [alignment.sh]

                quality [quality.sh]

            sample...1

                sample_analysis
		
		    qual [results of quality analysis]
		    
		    tmp [temporary files]
		    
		sample… 2

                    ...
---
