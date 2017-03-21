#!perl

###Don't use single end reads in SPAdes

while(my $log_file = glob "log.txt"){
	if($log_file eq 'log.txt'){
   open(INN, "log.txt") or die"";
   $/='>>';
   <INN>;
   while(<INN>){
   	chomp;
   	($name,$log_infor) = split(/\n/,$_,2);
   	$name =~ s/Assemble//g;
   	$name =~ s/\s+//g;
   	if($log_infor=~/Thank/){
   		$finishdb{$name} += 1;
   		}
   	}
   close IN;		
	}else{
		next;
		}
	}
	
###$ARGV[0] input library information, example format
#>HA1-4-CP10_CGTACTAG-CTCTCTAT
#Mix-1/HA1-4-CP10_CGTACTAG-CTCTCTAT_L001_R1.fastq.gz     /share/workplace/data/All_project/2_Sugarcane_BAC/Mix-1/HA1-4-CP10_CGTACTAG-CTCTCTAT_L001_R1.fastq.gz
#Mix-1/HA1-4-CP10_CGTACTAG-CTCTCTAT_L001_R2.fastq.gz     /share/workplace/data/All_project/2_Sugarcane_BAC/Mix-1/HA1-4-CP10_CGTACTAG-CTCTCTAT_L001_R2.fastq.gz

open(IN, $ARGV[0]) or die""; 
$/='>';
<IN>;
while(<IN>){
	chomp;
	($Lib, $lines) = split(/\s+/,$_,2);
	next if(exists($finishdb{$Lib}));   ####pass if the lib has been assembled before!
	print ">> Assemble $Lib\n";
	my $t1 = localtime(time());
  print "Starting time: $t1\n";
	mkdir $Lib;
	@data = split(/\n/,$lines);
################################
	print "1. Preparing data...\n";
	foreach $i(0..$#data){
		($reads, $path) = split(/\s+/,$data[$i]);
		if($reads =~ /(Mix.*)\/(.*.fastq.gz)/){
			$reads = $1."_".$2;
			}
		$cmd = "cp $path $Lib/$reads";
		system($cmd);
		}
	$cmd = "gunzip $Lib/*.gz"; 
	print "2. Uncompressing data ...\n"; 
	system($cmd);
################################	
	print "3. Combining data from different lanes ... \n";
	$cmd = "cat $Lib/*R1* > $Lib/reads.1.raw.fq"; system($cmd);
	$cmd = "cat $Lib/*R2* > $Lib/reads.2.raw.fq"; system($cmd);
	$cmd = "rm $Lib/*.fastq &";  system($cmd);
################################	
	print "4. Running Trimmomatic ...\n";
	$cmd  = "java -jar /home/zhangxt/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 10 -phred33 ";
	$cmd .= "$Lib/reads.1.raw.fq $Lib/reads.2.raw.fq ";
	$cmd .= "$Lib/reads.1.pe.fq $Lib/reads.1.se.fq ";
	$cmd .= "$Lib/reads.2.pe.fq $Lib/reads.2.se.fq ";
	$cmd .= "ILLUMINACLIP:/home/zhangxt/software/Trimmomatic-0.33/adapters/adapter.fasta:2:30:10 LEADING:3 TRAILING:3 ";
	$cmd .= "SLIDINGWINDOW:4:15 MINLEN:50";
	system($cmd);
################################	
	print "5. Removing duplication ...\n";
	open(OUT, ">$Lib/file.list") or die"";
	print OUT "$Lib/reads.1.pe.fq\n";
	print OUT "$Lib/reads.2.pe.fq\n";
	$cmd = "fastuniq -i $Lib/file.list -t q -o $Lib/reads.1.rdup.fq -p $Lib/reads.2.rdup.fq";
	system($cmd);
################################		
  $cmd = "cat $Lib/reads*.se.fq > $Lib/readsSE.fq";  system($cmd);
  $cmd = "rm $Lib/*.raw.fq "; system($cmd);
  $cmd = "rm $Lib/*.se.fq "; system($cmd);
  $cmd = "rm $Lib/*.pe.fq";  system($cmd);
  
################################	  
  print "6. Aligning reads to contamination...\n"; 
  $cmd = "cat $Lib/*rdup.fq > $Lib/fastq_file.fastq"; system($cmd);  
  $cmd = "parallel-fastq -j 12 \"lastal -Q1 -e120  contamination/condb \" < $Lib/fastq_file.fastq >  $Lib/myalns.maf ";
  system($cmd);
################################
  print "7. Removing contamination...\n"; 	 
	$cmd1 = "perl bin/remove_contamination.pl $Lib/myalns.maf $Lib/reads.1.rdup.fq > $Lib/reads.1.clean.fq ";
	$cmd2 = "perl bin/remove_contamination.pl $Lib/myalns.maf $Lib/reads.2.rdup.fq > $Lib/reads.2.clean.fq ";
	$cmd  = "perl bin/remove_contamination.pl $Lib/myalns.maf $Lib/readsSE.fq > $Lib/readsSE.clean.fq ";
	system($cmd1); system($cmd2); system($cmd);
	$cmd = "rm $Lib/readsSE.fq"; system($cmd);
	$cmd = "rm $Lib/fastq_file.fastq"; system($cmd);
	$cmd = "rm list.txt"; system($cmd);
  $cmd = "rm $Lib/myalns.maf"; system($cmd);
  $cmd = "rm $Lib/*rdup.fq"; system($cmd);
################################	  
  print "8. SPAdes assembly ...\n";
  $cmd = "spades.py -t 20 --careful -1 $Lib/reads.1.clean.fq -2 $Lib/reads.2.clean.fq -o $Lib/OUT";
  system($cmd);
################################	  
  print "9. check Assembly ...\n";
  $scaf_path = $Lib."/OUT/scaffolds.fasta";
  print "Statistics for Assembly results: \n";
  $cmd = "perl ~/software/script/faSize.pl $scaf_path "; system($cmd);
  my $t2 = localtime(time());
  print "Ending time: $t2\n";
	}
close IN;


