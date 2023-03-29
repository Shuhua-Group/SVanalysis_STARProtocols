#!/usr/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec;

#==========
#==========
my $help;

my $dir;
my $srr_list;
my $out_prefix;
my $input_suffix;


GetOptions(
    'dir|d=s' => \$dir,
    'srr_list|s=s' => \$srr_list,
    'out_prefix|o=s' => \$out_prefix,
    'input_suffix|i=s' => \$input_suffix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  IGVtools/code_bam2xml.pl

  Options:
  --dir|d <STR>
  --srr_list|s <STR>        
  --out_prefix|o <STR>
  --input_suffix|i <STR>
  --help or -h         output help message

=cut


# my @suffix_list=qw(.vcf.gz .vcf);
# my ($base, $dir, $suffix) = fileparse($vcf, @suffix_list);

my $basename_out=$out_prefix.".htseq.gene.count.mat";
my $outname = File::Spec->catfile($dir,$basename_out);
open(OUT,">$outname");

open(IN,$srr_list);
my %hash_name_arr;
my %hash_name_id;
my @sampleName=();
while(<IN>)
{
    chomp;
    my $basename=$_.".".$input_suffix;
    my $filename = File::Spec->catfile($dir,$basename);
    push(@sampleName,$_);
    open(IN1,$filename);
    ## all sample have the same gene number
    while(<IN1>)
    {
        chomp;
        if(/^\_\_/)
        {
            next;
        }
        my ($name,$id,$num)=split(/\s+/,$_);
        if(exists $hash_name_id{$name})
        {
            push(@{$hash_name_arr{$name}},$num);
        }else
        {
            $hash_name_id{$name}=$id;
            @{$hash_name_arr{$name}}=($num);
        }        
    }
    close IN1;
}


print OUT "Gene\tSymbol";
foreach my $i(@sampleName)
{
    print OUT "\t$i";
}
print OUT "\n";

foreach my $name (keys %hash_name_id)
{
    print OUT "$name\t$hash_name_id{$name}";
    foreach my $i(@{$hash_name_arr{$name}})
    {
        print OUT "\t$i";
    }
    print OUT "\n";

}
close OUT;

# sub save_in_array{
#     my ($hash_ref,$value)=@_;
#     if(@{$hash_ref}==0)
#     {
#         @{$hash_ref}=();
#     }
#     push(@{$hash_ref},$value);
#     # if($value eq "None")
#     # {
#     #     return;
#     # }
#     # my %tmp=();
#     # map{$tmp{$_}=1}@{$hash_ref};
#     # if(not exists $tmp{$value})
#     # {
#     #     push(@{$hash_ref},$value);
#     # }
# }

