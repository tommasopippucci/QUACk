##!/usr/bin/perl
#!/apps/perl/5.26.1

use warnings;
use strict;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(any uniq);
use Getopt::Long qw(GetOptions);
use Path::Tiny;

# Author: TP Aug 2019 
# Name: QUACk (QUality Alignment and Calling)
# Aim: To prepare bash scripts to run NGS pipeline 
# Usage:quack.pl -sample_sheet <SAMPLE_FILE.sample> -param <PARAMETER_FILE.param> -reference <hg19|hg38|..> -project <PATH/TO/PROJECT/FOLDER> [-aligner <default:BWA> -threads <default:1> -fastq_extension <default:fastq.gz>]

                                            #Declaring cycling $ variables
my($i) = 0; my($j) = 0; my($k) = 0; my $bam;
                                            #Declaring input argv $ variables
my($sample_file); my($param_file); my($threads)=1; my($project); my($ref); my($align)="bwa"; my($fastq_ext)="fastq.gz"; my($suffix)=""; my ($bed_ext)="";
                                            #Declaring sample sheet $ variables
my($sample); my($flow_cell); my($lane); my($index); my($enrichment); my($target_set); my($library); my($platform); my($provider);
                                            #Declaring parameter file $ variables
my($ref_file); my($target_dir); my($align_cmd); my($samtools_cmd); my($picard_cmd); my($gatk_cmd); my($gatk3_cmd); my($working_dir); my($fastq_dir); my($gvcf_dir); my($bam_dir); my($operator); my($id); my($refseq_genes);
                                            #Declaring date and time $ variables
my($second,$minute,$hour,$monthday,$month,$year,$weekday,$yearday) = localtime(time); my($datetime);
                                            #Declaring user name $ variables
my $username = $ENV{LOGNAME} || getpwuid($<) || $ENV{USER};
                                            #Declaring @ variables
my(@line) = (); my(@month_names)  = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec); my(@bam) = (); my(@known_sites) = ();
                                            #Declaring % variables
my(%sample_info);

                                                                                         #Expressing date and time
$year += 1900;
$year = sprintf("%02d", $year % 100);
$datetime = "$month_names[$month]$monthday$year-$hour.$minute.$second";

print "You ($username) are running QUACk (QUality Alignment and Calling) on $datetime\n";
print "Usage:quack.pl -sample_file <SAMPLE_FILE.sample> -parameter_file <PARAMETER_FILE.param> -reference <hg19|hg38|..> -project <PATH/TO/PROJECT/FOLDER> [-aligner <default:BWA> -threads <default:1> -fastq_extension <default:fastq.gz>]\n";

GetOptions( 'sample_file=s' => \$sample_file,
            'parameter_file=s' => \$param_file,
            'project=s' => \$project,
            'reference=s' => \$ref,
            'aligner=s' => \$align,          
            'threads=i' => \$threads,
            'fastq_extension=s' => \$fastq_ext,
            'fastq_suffix=s' => \$suffix,
            'extended_bed=s' => \$bed_ext
) or die "Invalid arguments!";

die "Missing sample sheet specification" unless $sample_file;
die "Missing parameter file specification" unless $param_file;
die "Missing reference sequence specification" unless $param_file;
die "Missing project name specification" unless $project;

                                                                                         #Checking and creating directories
open PARAM, "$param_file" or die $!;

while (<PARAM>)
   {
   next if ($_ =~ /^\#/);
   chomp;
#tolto da Fede&Marta
#s/\r//;s/\s+//;
   @line=split(",", $_);
   if ($line[0] eq $ref)
      {
      $ref_file=$line[1];
      }
   if ($line[0] eq "target_dir")
      {
      $target_dir=$line[1];
      }
   if ($line[0] eq "target_set")
      {
      $target_set=$line[1];
      }
   if ($line[0] eq $align)
      {
      $align_cmd=$line[1];
      }
   if ($line[0] eq "working_dir")
      {
      $working_dir=$line[1];
      }
   if ($line[0] eq "fastq_dir")
      {
      $fastq_dir=$line[1];
      }
   if ($line[0] eq "bam_dir")
      {
      $bam_dir=$line[1];
      }
   if ($line[0] eq "gvcf_dir")
      {
      $gvcf_dir=$line[1];
      }      
   if ($line[0] eq "known_sites")
      {
      @known_sites=@line[1..$#line];
      }
   if ($line[0] eq "samtools")
      {
      $samtools_cmd=$line[1];
      }
   if ($line[0] eq "picard")
      {
      $picard_cmd=$line[1];
      }
   if ($line[0] eq "gatk")
      {
      $gatk_cmd=$line[1];
      }
   if ($line[0] eq "gatk3")
      {
      $gatk3_cmd=$line[1];
      }
   if ($line[0] eq "refseq_genes")
      {
      $refseq_genes=$line[1];
      }
   }
   
close PARAM;

open SAMPLE, "$sample_file" or die $!;

while (<SAMPLE>)
   {
   next if ($_ =~ /^\#/);
   chomp;s/\r//;
   @line=split("\t", $_);
   if ($. == 1)
      {
      for $i (0..$#line)
         {
         if ($line[$i] eq "Sample") {$sample=$i}; if ($line[$i] eq "FlowCell") {$flow_cell=$i}; if ($line[$i] eq "Lane") {$lane=$i}; if ($line[$i] eq "Index") {$index=$i}; if ($line[$i] eq "Enrichment") {$enrichment=$i};
         if ($line[$i] eq "Library") {$library=$i}; if ($line[$i] eq "Platform") {$platform=$i}; if ($line[$i] eq "Provider") {$provider=$i}; 
         }
      }
   else
      {
      $id = $line[$sample]."_".$line[$flow_cell]."_L".$line[$lane];
      $sample_info{$line[$sample]}{$id} = [ @line ];
      }
   }

close SAMPLE;

if (-d "$working_dir/$project")
   {
   print "Project directory already exists. No need to create it\n";
   }
else
   {
   mkdir("$working_dir/$project") or die $!;
   print "Created project directory\n";
   }

if (-d "$bam_dir/$project")
   {
   print "BAM project directory already exists. No need to create it\n";
   }
else
   {
   mkdir("$bam_dir/$project") or die $!;
   print "Created BAM project directory\n";
   }

if (-d "$gvcf_dir/$project")
   {
   print "GVCF project directory already exists. No need to create it\n";
   }
else
   {
   mkdir("$gvcf_dir/$project") or die $!;
   print "Created GVCF project directory\n";
   }

mkdir("$working_dir/$project/$project\_$datetime") or die $!;
mkdir("$working_dir/$project/$project\_$datetime/bash/") or die $!;
mkdir("$working_dir/$project/$project\_$datetime/bash/quality/") or die $!;
mkdir("$working_dir/$project/$project\_$datetime/bash/aligncall/") or die $!;
mkdir("$working_dir/$project/$project\_$datetime/bash/coverage/") or die $!;

open LOG, "+>$working_dir/$project/$project\_$datetime/$project\_$datetime.log" or die $!;

print LOG "QUACk\n";
print LOG "User: $username\n";
print LOG "Date: $datetime\n";
print LOG "Command-line arguments:\n\tsample_file: $sample_file\n\tparameter_file: $param_file\n\tproject: $project\n\treference: $ref\n\t[aligner: $align, threads: $threads, fastq_extension: $fastq_ext]\n";
print LOG "Working directory: $working_dir\n";
print LOG "BAM files directory: $bam_dir\n";
print LOG "GVCF files directory: $gvcf_dir\n";
print LOG "Target file path: $target_dir/$target_set\n";

foreach $sample (keys %sample_info)
   {
   if (-d "$bam_dir/$project/$sample")
      {
      print "BAM project $sample directory already exists. No need to create it\n";
      }
   else
      {
      mkdir("$bam_dir/$project/$sample") or die $!;
      print "Created BAM project $sample directory\n";
      }
   if (-d "$gvcf_dir/$project/$sample")
      {
      print "GVCF project $sample directory already exists. No need to create it\n";
      }
   else
      {
      mkdir("$gvcf_dir/$project/$sample") or die $!;
      print "Created GVCF project $sample directory\n";
      }
   mkdir ("$working_dir/$project/$project\_$datetime/$sample") or die $!;
   mkdir ("$working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime") or die $!;
   mkdir ("$working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/") or die $!;
   mkdir ("$working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/") or die $!;

   open QUALITY, "+>$working_dir/$project/$project\_$datetime/bash/quality/pipeline_quality\_$sample.$datetime.sh" or die $!;
   open ALIGNCALL, "+>$working_dir/$project/$project\_$datetime/bash/aligncall/pipeline_aligncall\_$sample.$datetime.sh" or die $!;
   open COVERAGE, "+>$working_dir/$project/$project\_$datetime/bash/coverage/pipeline_coverage\_$sample.$datetime.sh" or die $!;
   
   print QUALITY  "#!/bin/bash\n\n";
   print ALIGNCALL "#!/bin/bash\n\n";
   print COVERAGE "#!/bin/bash\n\n";
   
   foreach $id ( sort keys %{ $sample_info{$sample} } )
      {
      print LOG "Enrichment kit used for $id: $sample_info{$sample}{$id}[$enrichment]\n";
      print QUALITY "#Performing fastQC analysis on fastQ files\n\n";
      print LOG "fastqc -v\n";
      print 
      print QUALITY "fastqc -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/ $fastq_dir/$sample/$id\_$sample_info{$sample}{$id}[$index]\_R1$suffix.$fastq_ext\n";
      print QUALITY "fastqc -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/ $fastq_dir/$sample/$id\_$sample_info{$sample}{$id}[$index]\_R2$suffix.$fastq_ext\n";
      
      print ALIGNCALL "#Creating BAM file $id\_$sample_info{$sample}{$id}[$index].bam, sorting and indexing BAM file $id\_$sample_info{$sample}{$id}[$index].sort.bam\n\n";           
      print ALIGNCALL "$align_cmd mem -M -t $threads -R \'\@RG\\tID:$sample\_$sample_info{$sample}{$id}[$lane]\_$sample_info{$sample}{$id}[$flow_cell]\\tSM:$sample\\tDT:$datetime\\tLB:$sample_info{$sample}{$id}[$library]\\tPL:$sample_info{$sample}{$id}[$platform]\\tCN:$sample_info{$sample}{$id}[$provider]\' $ref_file $fastq_dir/$sample/$id\_$sample_info{$sample}{$id}[$index]\_R1$suffix.$fastq_ext $fastq_dir/$sample/$id\_$sample_info{$sample}{$id}[$index]\_R2$suffix.$fastq_ext | $samtools_cmd sort -n -@ $threads -T $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/ -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].bam -\n";
      print ALIGNCALL "$samtools_cmd fixmate -m -@ $threads $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].bam $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].fixmate.bam\n";
     print ALIGNCALL "$samtools_cmd sort -@ $threads -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].fixmate.sort.bam $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].fixmate.bam\n"; 
      print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].fixmate.bam\n";
      print ALIGNCALL "$samtools_cmd index $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$id\_$sample_info{$sample}{$id}[$index].fixmate.sort.bam\n\n";
      $bam = $working_dir."/".$project."/".$project."_".$datetime."/".$sample."/".$sample."_".$datetime."/tmp/".$id."_".$sample_info{$sample}{$id}[$index].".fixmate.sort.bam";
      push @bam, $bam;
      }
   if (scalar @bam>1)
      {
      print ALIGNCALL "#Merging multiple BAM file into BAM file $sample\_$datetime/tmp/$sample.bam, sorting and indexing $sample.sort.bam\n\n";
      
      print ALIGNCALL "$samtools_cmd merge $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.bam";
      for $i (0..$#bam)
         {
         print ALIGNCALL " $bam[$i]";
         }
      print ALIGNCALL "\n";
      for $i (0..$#bam)
         {
         print ALIGNCALL "rm $bam[$i]\n";
         print ALIGNCALL "rm $bam[$i].bai\n";
         }
      print ALIGNCALL "$samtools_cmd sort -@ $threads -T $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/ -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.bam\n";
      print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.bam\n";
      print ALIGNCALL "$samtools_cmd index $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam\n\n";
      }
   elsif (scalar @bam == 1)
      {
      print ALIGNCALL "#Converting sorted BAM file name into $sample.sort.bam\n\n";
      
      print ALIGNCALL "mv $bam[0] $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam\n";
      print ALIGNCALL "$samtools_cmd index $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam\n\n";
      }
   else
      {
      die "No BAM file. Abort!\n"
      }

   print ALIGNCALL "#Marking duplicates in BAM file $sample.sort.bam\n\n";
   
   print ALIGNCALL "$samtools_cmd markdup-s -f $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/$sample.sort.markdup.txt -T $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/ -@ $threads $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam\n";
   print ALIGNCALL "$samtools_cmd index $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam\n\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.bam.bai\n";
   
   print ALIGNCALL "#Producing Samtools flag statistics for BAM file $sample.sort.markdup.bam\n\n";
   
   print ALIGNCALL "$samtools_cmd flagstat $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam > $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/$sample.sort.markdup.flagstat.txt\n\n";
   
   print ALIGNCALL "#Checking BAM file $sample.sort.bam\n\n";
   
   print ALIGNCALL "$samtools_cmd quickcheck -v $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam > $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/$sample.sort.markdup.quickcheck.txt\n\n";
   
   print ALIGNCALL "#Applying GATK Base Quality Score Recalibration on BAM file $sample.sort.markdup.bam, sorting and indexing BAM file $sample.sort.markdup.recal.bam\n\n";
   
   print ALIGNCALL "$gatk_cmd BaseRecalibrator -I $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam -O $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal_data.table -R $ref_file";
   for $i (0..$#known_sites)
      {
      print ALIGNCALL " --known-sites $known_sites[$i]";
      }
   print ALIGNCALL "\n";

   print ALIGNCALL "$gatk_cmd ApplyBQSR -R $ref_file -I $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam --bqsr-recal-file $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal_data.table -O $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal.bam\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bam\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.bai\n";
   print ALIGNCALL "$samtools_cmd sort -T $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/ -o $bam_dir/$project/$sample/$sample.markdup.recal.sort.bam $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal.bam\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal.bam\n";
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal.bai\n";
#modificato da Fede&Marta
   print ALIGNCALL "rm $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/tmp/$sample.sort.markdup.recal_data.table\n";
###############################
   print ALIGNCALL "$samtools_cmd index $bam_dir/$project/$sample/$sample.markdup.recal.sort.bam\n\n";

   print ALIGNCALL "#Generating gVCF $sample.g.VCF from BAM file $sample.markdup.recal.sort.bam with GATK\n\n";

   print ALIGNCALL "$gatk_cmd HaplotypeCaller -R $ref_file -I $bam_dir/$project/$sample/$sample.markdup.recal.sort.bam --emit-ref-confidence GVCF -L $target_dir/$target_set$bed_ext.bed -O $gvcf_dir/$project/$sample/$sample.g.vcf\n";
   
   print COVERAGE "$gatk3_cmd -T DepthOfCoverage -R $ref_file -I $bam_dir/$project/$sample/$sample.markdup.recal.sort.bam -L $target_dir/$target_set.bed -geneList $refseq_genes -o $working_dir/$project/$project\_$datetime/$sample/$sample\_$datetime/qual/$sample.coverage -omitBaseOutput -ct 20\n";
   
   @bam = ();

   close QUALITY;
   close ALIGNCALL;
   close COVERAGE;
   }
   
close LOG;
