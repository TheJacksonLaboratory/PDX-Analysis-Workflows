#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<'USAGE';

#########################GeneName_and_Normalization_without_UCSC.pl#############################

        usage: GeneName_and_Normalization_without_UCSC.pl [options]

               -i1 = infile1    [RSEM ENSEMBL genes.results]
               -i2 = infile2    [RSEM ENSEMBL isoforms.results]
               -a1 = accession1 [accession File1]

##################################################################################

USAGE

my ($infile1, $infile2, $accessionFile1);

my $result = GetOptions("i1=s"  => \$infile1, "i2=s" => \$infile2, "a1=s"  => \$accessionFile1);

die $usage unless ($infile1);            ##### Mandatory arguments
die $usage unless ($infile2);            ##### Mandatory arguments
die $usage unless ($accessionFile1);     ##### Mandatory arguments

###########################################################################################################################


my $seventyFifthQuartileGene  = quantileCal($infile1); 
my $seventyFifthQuartileIso   = quantileCal($infile2); 

#######Gene.results modification#############



open(FILEACC, $accessionFile1) || die "cannot open the $accessionFile1 ";####hg38/mm10_final

my %hashGene = ();

while(my $readFileACC = <FILEACC>)
{
  if($readFileACC =~ /^\s*.*?\s+(.*?)\s+(.*?)\s+(.*)$/)
  {
    my $key1 = $1;
    my $valueG = "$2";
#   my $UCSC = $3;

    if(!exists($hashGene{$key1}))
    {
     $hashGene{$key1}= $valueG; ####making ID by ENSG
    }

  }

} 


open(FILEINGENE, $infile1) || die "cannot open the $infile1 file";####*genes.results 
open(FILEOUTGENEName, ">$infile1.withGeneName") || die "cannot open the file"; 
open(FILEOUTGENENorm, ">$infile1.Normalized")   || die "cannot open the file";


my $flagGene = 0;
my ($keyG, $ENSName, $ExpectedcountG);

while(my $readFileGene = <FILEINGENE>)
{ 
  if($flagGene == 0 )
  {
    $flagGene = 1;
    chomp $readFileGene;
    print FILEOUTGENEName "$readFileGene\tGeneName\n";
    print FILEOUTGENENorm "gene_id\ttranscript_id(s)\tGeneName\tnormalized_count\n";
    next;
  }

  if($flagGene == 1)
  {
    if($readFileGene =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
    { 
       $keyG              = $1;
       $ENSName           = $2;
       $ExpectedcountG = sprintf("%.2f",(($5/$seventyFifthQuartileGene)*1000));
      
      if(exists($hashGene{$keyG}))
      {
        chomp $readFileGene;
        print FILEOUTGENEName "$readFileGene\t$hashGene{$keyG}\n";
        print FILEOUTGENENorm "$keyG\t$ENSName\t$hashGene{$keyG}\t$ExpectedcountG\n";
      } 

      else
      {
        chomp $readFileGene;
        print FILEOUTGENEName "$readFileGene\t-\n";
        print FILEOUTGENENorm "$keyG\t$ENSName\t-\t$ExpectedcountG\n";

      }


    }
  }
}


##################################################################################################################


close(FILEACC);

#######Isoforms.results modification#############


open(FILEACC, $accessionFile1) || die "cannot open the $accessionFile1 ";####hg19/mm10_final

my %hashIso = ();

while(my $readFileACC = <FILEACC>)
{
  if($readFileACC =~ /^\s*(.*?)\s+.*?\s+(.*?)\s+(.*)$/)
  {
    my $key2 = $1;
    my $valueI = "$2";
    #my $UCSC_I = $3;

    if(!exists($hashIso{$key2}))
    {
     $hashIso{$key2}= $valueI; ####making ID by ENST
    }

  }

} 


open(FILEINISO, $infile2) || die "cannot open the $infile2 file";####*genes.results 
open(FILEOUTISOName, ">$infile2.withGeneName") || die "cannot open the file"; 
open(FILEOUTISONorm, ">$infile2.Normalized")   || die "cannot open the file";

my $flagIso = 0;
my ($keyI, $ENSIName, $ExpectedcountI);

while(my $readFileIso = <FILEINISO>)
{ 
  if($flagIso == 0 )
  {
    $flagIso = 1;
    chomp $readFileIso;
    print FILEOUTISOName "$readFileIso\tGeneName\n";
    print FILEOUTISONorm "transcript_id\tgene_id\tGeneName\tnormalized_count\n";
    next;
  }

  if($flagIso == 1)
  {
    if($readFileIso =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
    { 
       $keyI               = $1;
       $ENSIName           = $2;
       $ExpectedcountI     = sprintf("%.2f",(($5/$seventyFifthQuartileIso)*300));
      
      if(exists($hashIso{$keyI}))
      {
        chomp $readFileIso;
        print FILEOUTISOName "$readFileIso\t$hashIso{$keyI}\n";
        print FILEOUTISONorm "$keyI\t$ENSIName\t$hashIso{$keyI}\t$ExpectedcountI\n";
      } 

      else
      {
        chomp $readFileIso;
        print FILEOUTISOName "$readFileIso\t-\n";
        print FILEOUTISONorm "$keyI\t$ENSIName\t-\t$ExpectedcountI\n";

      }


    }
  }
}



sub quantileCal
{

my $filename = shift;

open(FILEOUT1, ">$filename.TMP") || die "cannot open the file";

print FILEOUT1 "X=read.table(\"$filename\", header=T)\n";
print FILEOUT1 "head(X)\n";
print FILEOUT1 "Y=subset(X, X\$expected_count > 0)\n";
print FILEOUT1 "head(Y)\n";
print FILEOUT1 "Z=sort(Y[,5])\n";
print FILEOUT1 "quantile(Z)\n";


system("cat $filename.TMP | R --vanilla >$filename.TMP2 ");

open(FILEINTMP, "$filename.TMP2") || die "cannot open the file";

my $flag = 0;
my $value;

label:while(my $readFile = <FILEINTMP>)
{
 if($readFile =~ /quantile/)
  {
    $flag = 1; 
    next;
  }
  
 elsif($flag == 1)
  { 
    $flag = 2;
    next;
  }

 elsif($flag == 2)
  {
     if($readFile =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
     {
       $value = $4;
       last label;
     }

  }
}

close(FILEOUT1);
close(FILEINTMP);

system("rm $filename.TMP $filename.TMP2 ");

return $value;

}

