#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use lib "/public/home/PGG/liuqi/software";
use Soft; # use package by myself which include the path of software

#==========
my $config="/public/home/PGG/liuqi/software/config.txt";
my $bgzip=parse_config($config,"bgzip");
#==========
my $vcf='';
my $help;

GetOptions(
    'vcf|v=s' => \$vcf,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  code_convert_GT_FT.pl 

  Options:
  --vcf|v <STR>        graghtyper out vcf [mandatory]
  --help or -h         output help message
  
=cut

my $script_path = $Bin;
# 尖括号中间是搜索模式，尖括号运算符能返回与该模式匹配的文件列表，这称为一个glob，比如< *.bat>


my @suffix_list=qw(.vcf.gz .vcf);
my ($base, $dir, $suffix) = fileparse($vcf, @suffix_list);

my $out2=$base.".perSite.logfile";
open(OUT2,">$out2") or die "ERROR: $out2 can not be create!";
my $time=localtime();
print OUT2 "$time\tStart...\n";

#====================================
#  output individual's SV information
#====================================
my $out1=$base.".perSite.SVcount";
open(OUT1,">$out1") or die "ERROR: $out1 can not be create!";
print OUT1 "#chr\tstart\tend\tsv_name\tsv_type\tsv_size\tnum_of_ind_have_SV\tcount_of_sv_allele\tnum_of_ind_having_genotype\tnum_of_ind_inVCF\n";
#====================================
#  calculate individual's SV number
#====================================
if($suffix eq ".vcf.gz")
{
    open(IN,"gzip -dc $vcf | ") or die "ERROR: $vcf can not be open!";
}else
{
    open(IN,$vcf) or die $!;
}
while(<IN>)
{
    chomp;
    if(/^##/){next;}
    if(/^#C/)
    {
        next;
    }
    my @s=split(/\t/,$_);
    my $chr=$s[0];
    my $start=$s[1];
    my $name=$s[2];
    #  record basic information for SV
    my $alt=$s[4];
    my $info=$s[7]; 
    $info =~ /SVTYPE=(?<f1>[^;\s]+)/;
    my $sv_type=$+{f1};
    $info =~ /END=(?<f1>[^;\s]+)/;
    my $end=$+{f1};
    my $sv_size;
    if($info =~ /SVSIZE=(?<f1>[^;\s]+)/)
    {
        $sv_size=$+{f1};
    }elsif($alt =~ /SVSIZE=(?<f1>[^:\+]+)/)
    {
        $sv_size=$+{f1};
    }else
    {
        die "ERROR: $s[2] not have size";
    }
    my $ind_sum=$#s-8;
    print OUT1 "$chr\t$start\t$end\t$name\t$sv_type\t$sv_size";
    my $ind_num_1=0;
    my $count_1=0;
    my $ind_call=0;
    for(my $i=9;$i<=$#s;$i++)
    {
        my ($geno_i)=split(/:/,$s[$i],2);
        if($geno_i ne "./.")
        {
            $ind_call++;
            # this individual have genotype call
            if($geno_i =~ /1/)
            {
                $ind_num_1++;
                if($geno_i eq "1/1" || $geno_i eq "1|1")
                {
                    $count_1=$count_1+2;
                }else
                {
                    $count_1=$count_1+1;
                }
            }
        }
    }
    print OUT1 "\t$ind_num_1\t$count_1\t$ind_call\t$ind_sum\n";
}
close IN;
close OUT1;

$time=localtime();
print OUT2 "$time\tDone...\n";

