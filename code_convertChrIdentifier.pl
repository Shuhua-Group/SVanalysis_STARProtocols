#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

#==========
#==========
my $fasta;
my $relationship;
my $out_prefix;
my $help;

GetOptions(
    'fasta|f=s' => \$fasta,
    'relationship|r=s' => \$relationship,
    'out_prefix|o=s' => \$out_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  [Example] perl code_convertChrIdentifier.pl --fasta species.fasta --relationship relationship --out_prefix species.new.fasta

  Options:
  --fasta|f <STR>        reference genome fasta [mandatory]
  --relationship|r <STR> how to convert chromosome identifier
                  original chromosome ID must be the chromosome ID in reference genome fasta.
                  file format:
                  <original chromosome ID>\t<new chromosome ID>
                  e.g.
                  GWHAMMN00000001 1
  --out_prefix|o <STR>   prefix of out file
  --help or -h           output help message
  
=cut

open(IN,$relationship);
my %hash;
while(<IN>)
{
  chomp;
  my ($old,$new)=split(/\s+/,$_);
  $hash{$old}=$new;
}
close IN;

if($fasta =~ /gz$/)
{
    open(IN,"gzip -dc $fasta | ") or die $!;
}else
{
    open(IN,$fasta) or die $!;
}
my $out_file=$out_prefix.".fasta.gz";
open(OUT,"|gzip -c > $out_file") or die $!;

my $line=<IN>;
while($line)
{
  chomp($line);
  if($line =~ /^>/)
  {
    my ($fasta_name)=split(/\s+/,$line);
    $fasta_name=~s/>//;
    my $new_name=$hash{$fasta_name};
    print OUT ">$new_name\n";
    # print ">$fasta_name\t$new_name\n";
    $line=<IN>;
    while($line)
    {
      if($line =~ /^>/)
      {
        last;
      }else
      {
        print OUT "$line";
        $line=<IN>;
      }
    }
  }
}
close IN;
close OUT;

