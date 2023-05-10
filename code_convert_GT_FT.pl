#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

#==========
my $vcf='';
my $genotype_quality = 20;
my $depth=2;
my $help;

GetOptions(
    'vcf|v=s' => \$vcf,
    'gq|g=s' => \$genotype_quality,
    'dp|d=i' => \$depth,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  code_convert_GT_FT.pl 
  
  [Example] code_convert_GT_FT.pl --vcf input.vcf --gq 13 --dp 1
  
  Options:
  --vcf|v <STR>        graghtyper out vcf [mandatory]
  --gq|g <INT>        (>=gq)tgenotype_quality [default: 20]
  --dp|d <INT>        (>=dp)depth [defult: 2]
  --help or -h         output help message
  
=cut

my @suffix_list=qw(.vcf.gz .vcf);
my ($base, $dir, $suffix) = fileparse($vcf, @suffix_list);

if($suffix eq ".vcf.gz")
{
    open(IN,"gzip -dc $vcf | ") or die "$vcf can't be open!";
}else
{
    open(IN,$vcf) or die "$vcf can't be open!";
}
my $out_file=$base.".GQ".$genotype_quality.".DP".$depth.".vcf.gz";
open(OUT,"| $bgzip -c > $out_file") or die "$out_file can't be open!";
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AHRm.525
# (8)FORMAT:
# GT:FT:AD:MD:DP:RA:PP:GQ:PL
while(<IN>)
{
    chomp;
    if(/^#/)
    {
        print OUT "$_\n";
        next;
    }
    my @s=split(/\t/,$_);
    my $out1=join("\t",@s[0..8]);
    print OUT "$out1";
    for(my $i=9;$i<=$#s;$i++)
    {
        my @format=split(/:/,$s[$i]);
        my $dp=$format[4];
        my $gq=$format[7];
        if($dp<$depth || $gq<$genotype_quality)
        {
            $format[0]="./.";
        }
        # output format
        my $out2=join(":",@format);
        print OUT "\t$out2";
    }
    print OUT "\n";
}
close IN;
close OUT;
