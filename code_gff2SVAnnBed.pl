#!/bin/perl

#* This script is used to convert gff to a bed-format file, which is used in further annotation analysis.

#* [Example] perl code_gff2SVAnnBed.pl input.gff.gz out_prefix

my $gff=$ARGV[0];
my $prefix=$ARGV[1];

if($gff =~ /gz$/)
{
    open(IN,"gzip -dc $gff | ") or die $!;
}else
{
    open(IN,$gff) or die $!;
}
# bed
# [start,end)
my $out1=$prefix.".gff.cds.bed";
my $out2=$prefix.".gff.mRNA.bed";
my $out3=$prefix.".gff.UTR.bed";
my $out4=$prefix.".gff.promoter.bed";
my $out5=$prefix.".gff.upstream.bed";
my $out6=$prefix.".gff.downstream.bed";
my $out7=$prefix.".gff.gene.bed";

open(OUT1,">$out1") or die $!;
open(OUT2,">$out2") or die $!;
open(OUT3,">$out3") or die $!;
open(OUT4,">$out4") or die $!;
open(OUT5,">$out5") or die $!;
open(OUT6,">$out6") or die $!;
open(OUT7,">$out7") or die $!;

my $parent_id;
my $gene_name;
my %hash;
while(<IN>)
{
    chomp;
    if($_ =~ /^#/ or length($_)==0)
    {
        next;
    }
    my @s=split(/\t/,$_);
    my $attribute=$s[8];
    my $strand=$s[6];
    # 0-based bed
    # [start,end)
    $s[3]=$s[3]-1;
    if($s[2] eq "mRNA")
    {
        $attribute =~ /ID=(?<f1>[^;]+)/;
        $parent_id=$+{f1};
        $attribute =~ /Name=(?<f1>[^;]+)/;
        $gene_name=$+{f1};
        $hash{$parent_id}=$gene_name;    
        # mRNA
        # when SV not located in CDS but in mRNA, it will be intron
        print OUT2 "$s[0]\t$s[3]\t$s[4]\t${gene_name}:mRNA\tIntron\n";
        if($strand eq "+")
        {
            # + strand
            # promoter 1kb
            my $promoter_start=$s[3]-1000;
            if($promoter_start<0)
            {
                $promoter_start=0;
            }
            my $promoter_end=$s[3];
            print OUT4 "$s[0]\t$promoter_start\t$promoter_end\t${gene_name}:promoter\tPromoter\n";
            # upstream 20kb
            my $up20k_start=$s[3]-20000;
            if($up20k_start<0)
            {
                $up20k_start=0;
            }
            # my $up20k_end=$s[3];
            # 220618
            my $up20k_end=$s[3]-1000;
            print OUT5 "$s[0]\t$up20k_start\t$up20k_end\t${gene_name}:up20k\tUpstream20Kb\n";
            # upstream 50kb
            my $up50k_start=$s[3]-50000;
            my $up50k_end=$s[3]-20000;
            if($up50k_end>0)
            {
                if($up50k_start<0)
                {
                    $up50k_start=0;
                }
                print OUT5 "$s[0]\t$up50k_start\t$up50k_end\t${gene_name}:up50k\tUpstream50Kb\n";
            }
            # downstream 20kb
            my $down20k_start=$s[4];
            my $down20k_end=$down20k_start+20000;
            print OUT6 "$s[0]\t$down20k_start\t$down20k_end\t${gene_name}:down20k\tDownstream20Kb\n";
            # downstream 50kb
            my $down50k_start=$down20k_end;
            my $down50k_end=$down50k_start+50000;
            print OUT6 "$s[0]\t$down50k_start\t$down50k_end\t${gene_name}:down50k\tDownstream50Kb\n";
        }else
        {
            # - strand
            # promoter 1kb
            my $promoter_start=$s[4];
            my $promoter_end=$promoter_start+1000;
            print OUT4 "$s[0]\t$promoter_start\t$promoter_end\t${gene_name}:promoter\tPromoter\n";
            # upstream 20kb
            # my $up20k_start=$s[4];
            # 220618
            my $up20k_start=$promoter_end;
            my $up20k_end=$up20k_start+20000;
            print OUT5 "$s[0]\t$up20k_start\t$up20k_end\t${gene_name}:up20k\tUpstream20Kb\n";
            # upstream 50kb
            my $up50k_start=$up20k_end;
            my $up50k_end=$up50k_start+50000;
            print OUT5 "$s[0]\t$up50k_start\t$up50k_end\t${gene_name}:up50k\tUpstream50Kb\n";
            # downstream 20kb
            my $down20k_start=$s[3]-20000;
            if($down20k_start<0)
            {
                $down20k_start=0;
            }
            my $down20k_end=$s[3];
            print OUT6 "$s[0]\t$down20k_start\t$down20k_end\t${gene_name}:down20k\tDownstream20Kb\n";
            # downstream 50kb
            my $down50k_start=$s[3]-50000;
            my $down50k_end=$s[3]-20000;
            if($down50k_end>0)
            {
                if($down50k_start<0)
                {
                    $down50k_start=0;
                }
                print OUT6 "$s[0]\t$down50k_start\t$down50k_end\t${gene_name}:down50k\tDownstream50Kb\n";
            }   
        }
    }elsif($s[2] eq "CDS")
    {
        $attribute =~ /Parent=(?<f1>[^;]+)/;
        $parent_id2=$+{f1};
        $gene_name=$hash{$parent_id2};
        print OUT1 "$s[0]\t$s[3]\t$s[4]\t${gene_name}:CDS\tCDS\n";
    }elsif($s[2] =~ /UTR/)
    {
        $attribute =~ /Parent=(?<f1>[^;]+)/;
        $parent_id2=$+{f1};
        $gene_name=$hash{$parent_id2};
        print OUT3 "$s[0]\t$s[3]\t$s[4]\t${gene_name}:UTR\tUTR\n";
    }elsif($s[2] eq "gene")
    {
        $attribute =~ /Name=(?<f1>[^;]+)/;
        $gene_name=$+{f1};
        print OUT7 "$s[0]\t$s[3]\t$s[4]\t${gene_name}:gene\tGene\n";
    }else{}
}
close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
