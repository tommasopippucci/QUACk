# QUACk

QUACk (QUality Alignment and Calling) Perl script generating bash (.sh) scripts to perform three-steps NGS analysis (1. Read quality; 2. Alignment and alignment statistics; 3. Creation of gVCF)

## Usage

**Command**

    quack.pl 

**Mandatory options**

        -sample_file <SAMPLE_FILE.sample> 
        -parameter_file <PARAMETER_FILE.param> 
        -reference <hg19|hg38|..> 
        -project <PROJECT_FOLDER>
	
---
**Sample file** and **parameter file** are tab-delimited and comma-delimited files, and contain information about sample sequencing features and analysis parameters, respectively. 

**Reference** can be any released genome sequence pointing to a file path as specified in the **parameter file**. 

**Project folder** can be an existing or a non-existing directory where most analysis files will be located.

---

**Optional parameters**

        -aligner <default:bwa> 
        -threads <default:1> 
        -fastq_extension <default:fastq.gz> 
        -fastq_suffix <default:>
	
---

**Aligner** is the alignment tool name pointing to a command as specified in the **parameter file**. 

**Threads** is the number of processors to be allocated for tools allowing multi-threading.

**FastQ extension** is the term after the last dot in fastQ file name, the default term being *fastq.gz*.

**FastQ suffix** is an optional term that can be placed before the last dot in fastQ file name, the default being the *empty string*.

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


## Directory architectures

QUACk assumes that there are a *project* directory and its *sample* subdirectory(ies), named as the project argument in the command line and the sample name in the **sample file**, inside the *fastQ* directory as specified in the **parameter file**. Inside *sample* directories, QUACk will search for analysis input fastQ files.

---

    fastQ
        
	project
	
	    sample...1 [sample1 fastQ files]
	    
	    sample...2 [sample2 fastQ files]
	    
---

QUACk checks for the existence of the working, bam and gvcf directories.

Inside the *working* directory as specified in the **parameter file**, QUACk will generate a directory architecture where most final and all temporary analysis files will be located. First, QUACk will create a *project* directory. Inside this, every time it will be lauched QUACk will generate a new *project_analysis* directory flagged by the new analysis date and time. Inside every *project_analysis* directory there will be:
- a *log* file containing analysis parameter specifications
- a *bash* directory with *calling*, *alignment* and *quality* subdirectories each containing the *.sh* script pertaining to that part of the analysis
- *sample* directories each containing as many *sample_analysis* directories as the *project_analysis* directories are, inside which there will be the *tmp* directory for placing temporary files and the *qual* directory for storing different final files of the quality checks.

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

QUACk will also create a *project* directory in the *bam* and *gvcf* directories as specified in the **parameter file**, where respectively the final bam file and the gvcf file will be stored. While quality check and log files for every analysis will be mantained under the *working* directory architecture, final bam and gvcf files will be **OVERWRITTEN** at every new analysis.

---

    bam
        
	project
	
	    sample...1 [sample1 final bam file]
	    
	    sample...2 [sample1 final bam file]
	    
---

---

    gvcf
        
	project
	
	    sample...1 [sample1 gvcf file]
	    
	    sample...2 [sample1 gvcf file]
	    
---
