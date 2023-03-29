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
#==========
my $vcf = '';
my $min_sv_size = 50;
my $max_sv_size = 20000000;
my $help;

GetOptions(
    'vcf|v=s' => \$vcf,
    'len|l=i' => \$min_sv_size,
    'xlen|xl=i' => \$max_sv_size,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  code_SV_filtering_graphtyper2.pl

  Options:
   --vcf|v <STR>        svseq2 out vcf [mandatory]
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
my %hash_INV;
while(<IN>)
{
    chomp;
    if (/^#/)
    {
        next;
    }
    my @s=split(/\t/,$_);
    my $chr=$s[0];
    my $pos=$s[1];
    my $INFO=$s[7];
    $INFO =~ /SVTYPE=(?<f1>[^;\t]+)/;
    my $SVTYPE=$+{f1};
    $SVTYPE =~ s/[\<\>]//g;
    $INFO =~ /SVMODEL=(?<f1>[^;\t]+)/;
    my $SVMODEL=$+{f1};
    $INFO =~ /SVLEN=(?<f1>[^;\t]+)/;
    my $SVLEN=$+{f1};
    $SVLEN=~s/-//;
    if($SVTYPE eq "INV")
    {
        my $key=$chr.":".$pos.":".$SVLEN;
        $hash_INV{$key}{$SVMODEL}=1;
    }
}
close IN;

my $out=$base.".aggregated.vcf.gz";
open(OUT,"| bgzip -c > $out") or die "$out can't be open!";

if($suffix eq ".vcf.gz")
{
    open(IN,"gzip -dc $vcf | ") or die "$vcf can't be open!";
}else
{
    open(IN,$vcf) or die "$vcf can't be open!";
}
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
    $INFO =~ /SVMODEL=(?<f1>[^;\t]+)/;
    my $SVMODEL=$+{f1};
    if($SVTYPE eq "DEL" || $SVTYPE eq "DUP" || $SVTYPE eq "INS") 
    {
        if($SVMODEL eq "AGGREGATED")
        {
            if($SVTYPE eq "DEL" || $SVTYPE eq "DUP") 
            {
                $INFO =~ /SVLEN=(?<f1>[^;\t]+)/;
                my $SVLEN=$+{f1};
                $SVLEN=~s/-//;
                if($SVLEN >= $min_sv_size && $SVLEN <= $max_sv_size)
                {
                    print OUT "$_\n";
                }
            }elsif($SVTYPE eq "INS")
            {
                print OUT "$_\n";
            }
        }
    }elsif($SVTYPE eq "INV")
    {
        $INFO =~ /SVLEN=(?<f1>[^;\t]+)/;
        my $SVLEN=$+{f1};
        $SVLEN=~s/-//;
        if($SVLEN >= $min_sv_size && $SVLEN <= $max_sv_size)
        {
            my $chr=$s[0];
            my $pos=$s[1];
            my $key=$chr.":".$pos.":".$SVLEN;
            my $model="AGGREGATED";
            if($SVMODEL eq "AGGREGATED" && exists $hash_INV{$key}{$model})
            {
                # if INV have AGGREGATED for BREAKPOINT1 & BREAKPOINT2
                # output AGGREGATED
                print OUT "$_\n";
            }elsif(not exists $hash_INV{$key}{$model} && $SVMODEL ne "AGGREGATED")
            {
                # if INV only have BREAKPOINT1 or BREAKPOINT2
                # output BREAKPOINT1 or BREAKPOINT2
                print OUT "$_\n";
            }
        }
    }
}
close IN;
close OUT;
