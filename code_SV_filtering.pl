#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib $Bin;

#==========
my $tool = '';
my $vcf='';
my $rss = 3;
my $num_svs=2;
my $min_sv_size = 50;
my $max_sv_size = 20000000;
my $help;


GetOptions(
    'tool|t=s' => \$tool,
    'vcf|v=s' => \$vcf,
    'rss|r=s' => \$rss,
    'num|n=i' => \$num_svs,
    'len|l=i' => \$min_sv_size,
    'xlen|xl=i' => \$max_sv_size,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  code_SV_filtering.pl 

  Options:
  --tool|t <STR>       a tool name [mandatory]
  --vcf|v <STR>        svseq2 out vcf [mandatory]
  --rss|r <INT>        (>=RSS)the number of reads supporting the called SV allele [default: 3]
  --num|n <INT>        num of tools for joined_svs [defult: 2]
  --len|l <INT>        minimum size (bp) of DEL/DUP/INV to be kept [default: 50]
  --xlen|xl <INT>      maximum size (bp) of SV to be kept [default: 20000000]
  --help or -h         output help message

  Note: 
  tools can be: Manta, Lumpy,SVseq2,svimmer,graphtyper2
  For SVseq2, we retain DEL.
  For Manta, we retain DEL/DUP/INS/INV.
  For Lumpy, we retain DEL/DUP/INV.
  For svimmer, per sample, join Manta variants to other tools variants using svimmer and keep only variants that were also found by other tools
  For graphtyper2, we just keep <DEL:SVSIZE=XXX:AGGREGATED> following the paper's suggestion
  
=cut

my $script_path = $Bin;
# 尖括号中间是搜索模式，尖括号运算符能返回与该模式匹配的文件列表，这称为一个glob，比如< *.bat>
my @convert_script = <$script_path/code_SV_filtering_$tool*.pl>;
@convert_script = <$script_path/vcf_filter/code_SV_filtering_$tool*.pl> if (@convert_script == 0);


if (@convert_script == 0){
    die "TOOL name specified is absent in tools we analyzed or does not match the names:\n";
}
elsif (@convert_script > 1){
    die "TOOL name specified is redundant. Specify more strict name:\n";
}

# SV discovery
`$convert_script[0] --vcf $vcf --rss $rss --len $min_sv_size --xlen $max_sv_size` if ($tool eq 'Lumpy');
`$convert_script[0] --vcf $vcf --len $min_sv_size --xlen $max_sv_size` if ($tool eq 'Manta');
`$convert_script[0] --vcf $vcf --rss $rss --len $min_sv_size --xlen $max_sv_size` if ($tool eq 'SVseq2');
# SV merging
`$convert_script[0] --vcf $vcf --num $num_svs --len $min_sv_size --xlen $max_sv_size` if ($tool eq 'svimmer');

# SV graphtyper2 filtering
`$convert_script[0] --vcf $vcf --len $min_sv_size --xlen $max_sv_size` if ($tool eq 'graphtyper2');
