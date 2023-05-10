
#===============================================================================
#* This pipeline is used to genotype SV in all sample using graphtyper2
#TODO 140: SV genotyping in population-scale
#TODO 150: SV filtering followed by the recommendation by Graphtyper2
#===============================================================================

function usage()
{
    echo "
******This pipeline is used to genotype SV in all sample using graphtyper2
Usage: SVprotocol_part2.sh   <The long option is recommended>
    -s, --samplelist   (required) sample list, one sample per line
    -c, --config       (required) a config file including basic setting.
    -r, --repeat       (required) repeat sequence position.
    -j, --thread       (optional) the number of thread, default = 5
    -h, --help         (optional) This small usage guide

    1) Please prepare a config file following the example file
        1. Containing the path of software to be used
        2. You must use the variable name shown in the example file because, in this pipeline, we will use these variable names
    2) Please prepare a file including the repeat sequence region.
        0-based bed file.<chr><start><end><ID>
    3) Example:
    SVprotocol_140_150.sh --samplelist sample.list
                            --config example.config
                            --repeat example.repeat.bed
                            --thread 2
    "
}


OPTIONS=hs:c:j:r:f
LONGOPTIONS=help,samplelist:,config:,thread:,repeat:,force
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [ $? != 0 ]; then
    echo "ERROR: Terminating..."
    exit 2
fi
eval set -- "${PARSED}"

FLAG_sample=0
FLAG_config=0
FLAG_repeat=0
thread=5

while true
do
    case "$1" in
        -s|--samplelist)
            samplelist=$2
            FLAG_sample=1
            shift 2
            ;;
        -c|--config)
            config=$2
            FLAG_config=1
            shift 2
            ;;
        -r|--repeat)
            repeat=$2
            FLAG_repeat=1
            shift 2
            ;;
        -j|--thread)
            thread=$2
            shift 2
            ;;
        -h|--help)
            usage
            shift
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "ERROR: Unknown argument $1 . Please check again."
            exit 3
            ;;
    esac
done

until [ $FLAG_sample -eq 1 -a $FLAG_config -eq 1 ]
do
    echo "ERROR: missing required variables. Please check again."
    exit 1
done

function check_existed_file()
{
    local file_name=$1
    if [ -f $file_name ]
    then
        file_size=`ls -l $file_name |awk '{print $5}'`
        if [ $file_size -gt 0 ]
        then
            return 0
        else
            return 1
        fi
    else
        return 1
    fi
}

#===============================================================================
#****/ Main 
#===============================================================================

source $config
cd $your_analysis_dir

echo "SVprotocol_part2.sh"
printf "\t--samplelist\t%s\n" "$samplelist"
printf "\t--config\t%s\n" "$config"
printf "\t--thread\t%s\n" "$thread"
printf "%-20s\t%-s\n" "Hostname:" "`hostname`"
printf "%-20s\t%-s\n" "Working directory:" "`pwd`"
printf "INFO\t%-s\t%-s\n" "`date`" "Start"

#****/ Check required files
if [ ! -f ${ref_genome_fasta}.fai ]
then
    echo "Please index your fasta using bwa with command: bwa index ref.fasta"
	exit 1
fi


#************/ Create directory /********************************
mkdir -p 130.sv.detection/svimmer_Nsample >/dev/null 2>&1
mkdir -p 140.sv.genotyping >/dev/null 2>&1

#************/ Merging SVs discovered in all samples using svimmer software
svimmer_analysis_path=$your_analysis_dir/130.sv.detection/svimmer_Nsample
cd $svimmer_analysis_path
rm -f vcf.list >/dev/null 2>&1
cat $samplelist | while read sample
do
ln -s $your_analysis_dir/130.sv.detection/svimmer_1sample/$sample/${sample}.chr.svimmer.filtered.vcf.gz* ./ >/dev/null 2>&1
echo "${sample}.chr.svimmer.filtered.vcf.gz" >> vcf.list
done
#**** We extract chromosome name from genome.fasta.fai
#**** GWHAMMN00000001	325144201	17	60	61
chr=""
genome_fai=${ref_genome_fasta}.fai
while read line
do
    arr=($line)
    chr_i=${arr[0]}
    chr="$chr $chr_i"
done < $genome_fai
$python3_9 $svimmer --threads $thread vcf.list $chr | $bgzip -c > population.svimmer.vcf.gz
#**** Filtering SVs with sizes larger than 50bp and smaller than 2Mb
$perl ${this_protocol_script_path}/code_SV_filtering.pl --tool svimmer --vcf population.svimmer.vcf.gz --len 50 --xlen 2000000
$tabix -p vcf -f population.svimmer.filtered.vcf.gz

printf "INFO\t%-s\t%-s\n" "`date`" "svimmer_Nsample finished"

#************/ Subsampling regions with abnormally high sequence depth
graphtyper_analysis_path=$your_analysis_dir/140.sv.genotyping
cd $graphtyper_analysis_path
rm -f bam.list input.avg_cov_by_readlen >/dev/null 2>&1
cat $samplelist | while read sample
do
graphtyper_input_bam=${your_analysis_dir}/120.remove.duplicates/$sample/${sample}.dedup.bam
$samtools idxstats $graphtyper_input_bam | head -n -1 | awk '{sum+=$3+$4;ref+=$2}END{print sum/ref}' >> input.avg_cov_by_readlen
echo "$graphtyper_input_bam" >> bam.list
done
printf "INFO\t%-s\t%-s\n" "`date`" "Subsampling finished"

#************/	SV genotyping for each chromosome (5h / 12 CPUs)
#**** We extract chromosome name from genome.fasta.fai
#**** GWHAMMN00000001	325144201	17	60	61
genome_fai=${ref_genome_fasta}.fai
pop_svimmer=$svimmer_analysis_path/population.svimmer.filtered.vcf.gz
rm -f $graphtyper_analysis_path/vcf.list2 2>/dev/null
rm -f *aggregated.vcf.gz* 2>/dev/null
cat $genome_fai | while read line
do
arr=($line)
chr=${arr[0]}
start=1
end=${arr[1]}
cd $graphtyper_analysis_path
$graphtyper2 genotype_sv $ref_genome_fasta $pop_svimmer --sams=$graphtyper_analysis_path/bam.list --region=$chr:$start-$end --threads=$thread --force_use_input_ref_for_cram_reading --no_cleanup --max_files_open=200 --avg_cov_by_readlen=$graphtyper_analysis_path/input.avg_cov_by_readlen
#* Concatenating the result of graphtyper
cd $graphtyper_analysis_path/sv_results/$chr
vcf_list=vcf.list
rm -f $vcf_list 2>/dev/null
for vcf in `ls *vcf.gz`
do
    echo "$vcf" >> $vcf_list
done
out_vcf=$graphtyper_analysis_path/$chr.graphtyper2.vcf.gz
$bcftools concat --allow-overlap --file-list $vcf_list --output $out_vcf --output-type z --threads $thread
cd $graphtyper_analysis_path
#**** Filtering SVs with sizes larger than 50bp and smaller than 2Mb
$perl ${this_protocol_script_path}/code_SV_filtering.pl --tool graphtyper2 --vcf $out_vcf --len 50 --xlen 2000000 && $tabix -p vcf -f $chr.graphtyper2.aggregated.vcf.gz && rm -f $out_vcf >/dev/null 2>&1
echo "$graphtyper_analysis_path/$chr.graphtyper2.aggregated.vcf.gz" >> $graphtyper_analysis_path/vcf.list2
done
printf "INFO\t%-s\t%-s\n" "`date`" "graphtyper2 single chromosome:"

wait
#************/	Concatenating SV genotyping results in all chromosomes
cd $graphtyper_analysis_path
vcf_list=$graphtyper_analysis_path/vcf.list2
out_vcf=$graphtyper_analysis_path/population.graphtyper2.aggregated.vcf.gz
$bcftools concat --allow-overlap --file-list $vcf_list --output $out_vcf --output-type z --threads $thread && $tabix -p vcf -f $out_vcf

#**** remove temporary files
rm -rf $graphtyper_analysis_path/sv_results >/dev/null 2>&1
cat $vcf_list | xargs rm

printf "INFO\t%-s\t%-s\n" "`date`" "graphtyper2 all chromosome"

#************/	Filtering SVs followed by the recommendation by Graphtyper2
cd $graphtyper_analysis_path
filter_input_vcf=population.graphtyper2.aggregated
$vcffilter -f "( SVTYPE = DEL & QD > 12 & ( ABHet > 0.30 | ABHet < 0 ) & ( AC / NUM_MERGED_SVS ) < 25 ) | ( SVTYPE = DUP & QD > 5 & ( AC / NUM_MERGED_SVS ) < 25 ) | ( SVTYPE = INS & ( AC / NUM_MERGED_SVS ) < 25 & ( ABHet > 0.25 | ABHet < 0 ) & MaxAAS > 4 ) | ( SVTYPE = INV & ( AC / NUM_MERGED_SVS ) < 25  & ( ABHet > 0.25 | ABHet < 0 ) & MaxAAS > 4 )" ${filter_input_vcf}.vcf.gz | $bgzip -c > ${filter_input_vcf}.PASS.vcf.gz
#**** Converting genotype with DP<1 or GQ<13 as missing genotype
$perl ${this_protocol_script_path}/code_convert_GT_FT.pl --vcf ${filter_input_vcf}.PASS.vcf.gz --gq 13 --dp 1
$tabix -p vcf -f ${filter_input_vcf}.PASS.GQ13.DP1.vcf.gz

printf "INFO\t%-s\t%-s\n" "`date`" "graph filter finished"

#************/	Filtering SVs with missing rate<50%
graphtyper_out_pass_vcf=${filter_input_vcf}.PASS.GQ13.DP1.vcf.gz
filter_missing_rate_out=${graphtyper_out_pass_vcf/.vcf.gz/}.GENO0.5.vcf.gz
$vcftools --gzvcf $graphtyper_out_pass_vcf --max-missing 0.5 --recode --recode-INFO-all --stdout | $bgzip -c > $filter_missing_rate_out && $tabix -p vcf -f $filter_missing_rate_out

$perl ${this_protocol_script_path}/code_countINDnumPerSite.pl --vcf $filter_missing_rate_out

printf "INFO\t%-s\t%-s\n" "`date`" "missingRateL0.5 finished"

#************/	Filtering SVs overlapped with repeat region
SV_count=${filter_input_vcf}.PASS.GQ13.DP1.GENO0.5.perSite.SVcount
SV_bed=${SV_count/perSite.SVcount/bed}
#* SV_bed: 0-based
perl -alne 'if(/^#/){next;}$F[1]=$F[1]-1;print "$F[0]\t$F[1]\t$F[2]\t$F[0]-$F[1]-$F[2]-$F[5]-$F[4]-$F[3]"' $SV_count > $SV_bed
#* 2-overlap region
bedout_bed=population.repeat.bedout
toberemove_bed=population.repeat.rm.sv.bed
toberemove1=population.repeat.rm.list
# -wb: report region in A which overlaped with B, and original B
$bedtools intersect -a $SV_bed -b $repeat -wb > $bedout_bed && $perl ${this_protocol_script_path}/code_DUPrepeatRatio.pl $bedout_bed > $toberemove_bed && $perl -alne 'if($F[3]>=0.5){@s=split(/-/,$F[0]);print $s[$#s]}' $toberemove_bed > $toberemove1
#* 3-remove repeat SV
filter_repeat_out=${filter_missing_rate_out/.vcf.gz/}.rmRepeat.vcf.gz
$vcftools --gzvcf $filter_missing_rate_out --exclude $toberemove1 --recode --recode-INFO-all --stdout | $bgzip -c > $filter_repeat_out && $tabix -p vcf -f $filter_repeat_out
$perl ${this_protocol_script_path}/code_countINDnumPerSite.pl --vcf $filter_repeat_out


#************/ Ending
printf "INFO\t%-s\t%-s\n" "`date`" "Congratulations"
printf "INFO\t%-s\t%-s\n" "`date`" "End"


