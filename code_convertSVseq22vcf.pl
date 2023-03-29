#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
use Pod::Usage;


our vars($svseq2out_list,$out_prefix,$sample_name,$help);
GetOptions(
    "svseq2out_list|i=s"=>\$svseq2out_list,
    "out_prefix|o=s"=>\$out_prefix,
    "sample_name|s=s"=>\$sample_name,
    "help"=>\$help,
)
or pod2usage(-verbose => 0);
pod2usage(-verbose => 1) if (defined $help);
pod2usage(-verbose => 0) if ((! defined $svseq2out_list) and (! defined $out_prefix) and (! defined $sample_name));

=head1 SYNOPSIS

  code_convertSVseq22vcf.pl

  [Example] code_convertSVseq22vcf.pl -i svseq2.DEL.out.list -o svseq2

  Options:
   --svseq2out_list|i <STR>  <required> a list include svseq2 out files containing all chromosome
   --out_prefix|o <STR>      <required> output file's prefix
   --sample_name|s <STR>     <required> sample's name
   --help or -h             output help message
   
=cut



my $out4=$out_prefix.".logfile";
open(OUT4,">$out4") or die $!;

my $time=localtime();
print OUT4 "$time\tStart...\n";

# code come from https://github.com/stat-lab/EvalSVcallers

my %vcf;
my $flag = 0; # indicator which record the header

open(IN1,$svseq2out_list) or die $!;
while(<IN1>)
{
    chomp;
    if (!-e $_)
    {
        print STDERR "$_ is not found\n";
        next;
    }
    open (FILE, $_);
    my $type = '';
    my $pre_line = ''; # The information for chromosome, breakpoint, and SV size (breakpoint distance) was obtained from the space-separated header line for each called SV.
    my $reads = 0;
    while (my $line = <FILE>)
    {
        chomp $line;
        if ($line =~ /^#+$/)
        {
            if ($pre_line ne '')
            {
                # record previous SV information
                my @line = split (/\s+/, $pre_line);
                # range GWHAMMN00000010 477995 478010 478135 478150 0.937667 0.97
                # chr=GWHAMMN00000010
                # pos=477995
                # len=478135 - 477995 + 1
                my $len = 0;
                my $pos = 0;
                my $chr = '';
                my $end = 0;
                if ($type eq 'DEL')
                {
                    $chr = $line[1];
                    $pos = $line[2];
                    $len = $line[4] - $pos + 1;
                    $end = $line[4];
                }
                elsif ($type eq 'INS')
                {
                    $chr = $line[0];
                    $pos = $line[1];
                }
                my $chr02d = $chr;
                #$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
                ${${$vcf{$chr02d}}{$pos}}{$type} = "$len=$reads=$end";
            }
            $pre_line = '';
            $type = '';
            $flag = 0;
            $reads = 0;
        }
        elsif ($flag == 0)
        {
            # this line is header
            $pre_line = $line;
            $flag = 1;
            if ($line =~ /^range/)
            {
                # when line start with range, it mean that this is deletion
                $type = 'DEL';
            }
            else
            {
                # otherwise, this is insertion
                $type = 'INS';
            }
        }
        elsif ($flag == 1)
        {
            # other line is the read
            if ($type eq 'DEL')
            {
                $reads ++;
            }
            elsif ($type eq 'INS')
            {  
                # \d match number
                $reads ++ if ($line =~ /^\d+/);
            }
        }
    }
    close (FILE);
}
close IN1;

my $out_file=$out_prefix.".SVseq2.vcf";
open(OUT,">$out_file") or die $!;

print OUT "##fileformat=VCFv4.2\n";
print OUT "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print OUT "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print OUT "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
print OUT "##INFO=<ID=SU,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">\n";
print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";
my $id=1;
foreach my $chr (sort keys %vcf)
{
    foreach my $pos (sort {$a <=> $b} keys %{$vcf{$chr}})
    {
        foreach my $type (keys %{${$vcf{$chr}}{$pos}})
        {
            my ($len, $reads,$end) = split (/=/, ${${$vcf{$chr}}{$pos}}{$type});
            my $chr2 = $chr;
            #$chr2 =~ s/^0*//;
            print OUT "$chr2\t$pos\t$id\tN\t<$type>\t.\t.\tSVTYPE=$type;SVLEN=$len;SU=$reads;END=$end\tGT\t./.\n";
            $id=$id+1;
        }
    }
}
close OUT;

$time=localtime();
print OUT4 "$time\tend...\n";

close OUT4;
