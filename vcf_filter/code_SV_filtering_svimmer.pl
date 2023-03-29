#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
# use lib "/public/home/PGG/liuqi/software";
# use Soft; # use package by myself which include the path of software

#==========
# my $config="/public/home/PGG/liuqi/software/config.txt";
# my $bgzip=parse_config($config,"bgzip");
# #==========
my $vcf = '';
my $num_svs=2;
my $min_sv_size = 50;
my $max_sv_size = 20000000;
my $help;

GetOptions(
    'vcf|v=s' => \$vcf,
    'num|n=i' => \$num_svs,
    'len|l=i' => \$min_sv_size,
    'xlen|xl=i' => \$max_sv_size,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  code_SV_filtering_svimmer.pl

  Options:
   --vcf|v <STR>        svseq2 out vcf [mandatory]
   --num|n <INT>        num_svs [defult: 2]
   --len|l <INT>        minimum size (bp) of DEL/DUP/INV to be kept [default: 50]
   --xlen|xl <INT>      maximum size (bp) of SV to be kept [default: 20000000]
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

my $out=$base.".filtered.vcf.gz";
open(OUT,"| bgzip -c > $out") or die "$out can't be open!";
while(<IN>)
{
    chomp;
    if (/^#/)
    {
        print OUT "$_\n";
        next;
    }
    my @s=split(/\t/,$_);
    my $INFO=$s[7];
    $INFO =~ /SVTYPE=(?<f1>[^;\t]+)/;
    my $SVTYPE=$+{f1};
    $SVTYPE =~ s/[\<\>]//g;
    if($SVTYPE eq "DEL" || $SVTYPE eq "DUP" || $SVTYPE eq "INV") 
    {
        # num_joined_svs >=2
        $INFO =~ /SVLEN=(?<f1>[^;\t]+)/;
        my $SVLEN=$+{f1};
        $SVLEN=~s/-//;
        my $NUM_SVS=0;
        if($INFO =~ /NUM_JOINED_SVS=(?<f1>[^;\t]+)/)
        {
            $NUM_SVS=$+{f1};
        }elsif($INFO =~ /NUM_MERGED_SVS=(?<f1>[^;\t]+)/)
        {
            $NUM_SVS=$+{f1};
        }else
        {
            die "can't find NUM_*_SVS in $vcf\n";
        }
        if($SVLEN >= $min_sv_size && $SVLEN <= $max_sv_size)
        {
            if($NUM_SVS >=$num_svs)
            {
                print OUT "$_\n";
            }
        }
    }elsif($SVTYPE eq "INS")
    {
        # here, we only use INS from Manta, so num_joined_svs==1
        print OUT "$_\n";
    }
}
close IN;
close OUT;
