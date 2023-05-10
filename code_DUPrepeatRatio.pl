#!/usr/bin/perl


#* This script is used to calculate the cumulative length of region overlapped with repeats (1-bp overlapped threshold)
#* if the ratio of overlapped length > N, delete SV

#* [Example] perl code_DUPrepeatRatio.pl bedtools.bed.out

my $bedout=$ARGV[0];
# GWHAMMK00000148 82910   82931   GWHAMMK00000148-72601-122800-50199-DUP  GWHAMMK00000148 82910   82931   GWHAMMK00000148:82910:82931
open(IN,$bedout);
my $old=0;
my @old_seq=();
my $old_start=0;
my $old_len=0;
while(<IN>)
{
    chomp;
    my @s=split(/\t/,$_);  
    if($old ne $s[3])
    {
        # print old seq
        if($#old_seq>0)
        {
            my %d;
            $d{$_}++ for @old_seq;
            my $repeat_len=$d{1};
            my $repeat_ratio=$repeat_len/$old_len;
            print "$old\t$repeat_len\t$old_len\t$repeat_ratio\n";
        }
        # new seq
        
        my @raw_seq=split(/-/,$s[3]);
        $old=$s[3];
        $old_len=$raw_seq[3]+1;
        $old_start=$raw_seq[1];
        @old_seq=();

        # record first repeat
        $old_seq[0]=0;
        my $max_index=$raw_seq[3];
        $old_seq[$max_index]=0;

        my $bed_start=$s[1]-$old_start;
        my $bed_end=$s[2]-$old_start;
        my $bed_len=$bed_end-$bed_start+1;
        my @tmp=(1)x$bed_len;
        @old_seq[$bed_start..$bed_end]=@tmp;

        
    }else
    {
        my $bed_start=$s[1]-$old_start;
        my $bed_end=$s[2]-$old_start;
        my $bed_len=$bed_end-$bed_start+1;
        my @tmp=(1)x$bed_len;
        @old_seq[$bed_start..$bed_end]=@tmp;
    }
}
close IN;

# last record
# print old seq
if($#old_seq>0)
{
    my %d;
    $d{$_}++ for @old_seq;
    my $repeat_len=$d{1};
    my $repeat_ratio=$repeat_len/$old_len;
    print "$old\t$repeat_len\t$old_len\t$repeat_ratio\n";
}
