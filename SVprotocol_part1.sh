#!/usr/bin/sh

#===============================================================================
#* This pipeline is used to detect SV in one sample using three software
#TODO 110: mapping short-read sequence to reference genome
#TODO 120: Remove duplicates
#TODO 130: SV detection for one sample using three software
#===============================================================================

function usage()
{
    echo "
******This pipeline is used to detect SV in one sample using three software
Usage: SVprotocol_part1.sh   <The long option is recommended>
    -s, --sample   (required) sample id
    -c, --config   (required) a config file including basic setting.
    -j, --thread   (optional) the number of thread, default = 5
    -f, --force    (optional) removes output subdirectory if it exists, run all process again
    -h, --help     (optional) This small usage guide

    1) In your_analysis_dir, you must have a subdirectory named as fastq
        <sample_name> including two FASTQ file (<sample_name>_R1.fastq.gz and <sample_name>_R2.fastq.gz)
    2) Please prepare a config file following the example file
        1. Containing the path of software to be used
        2. You must use the variable name shown in the example file because, in this pipeline, we will use these variable names
    3) Please just include chromosome-level sequence or sequence with size larger than 250kb in the reference genome fasta file
    4) Before running this pipeline, please create dictionary for genome fasta
        >java -jar picard.jar CreateSequenceDictionary R=Rhipicephalus_microplus.chromosome.fasta O=Rhipicephalus_microplus.chromosome.dict
    5) ref_genome_fasta = <prefix>.fasta. Before you run SVseq2, please
        split genome-wide sequence fasta to single chromosome sequence fasta named <prefix>.<chr>.fasta, with index file named <prefix>.<chr>.fasta.fai
    6) Before running this pipeline, please index your reference genome by yourself
        >bwa index Rhipicephalus_microplus.chromosome.fasta
        Real time: 2288.064 sec; CPU: 2281.825 secs
        Rhipicephalus_microplus.chromosome.fasta
        Rhipicephalus_microplus.chromosome.fasta.amb
        Rhipicephalus_microplus.chromosome.fasta.ann
        Rhipicephalus_microplus.chromosome.fasta.bwt
        Rhipicephalus_microplus.chromosome.fasta.fai
        Rhipicephalus_microplus.chromosome.fasta.pac
        Rhipicephalus_microplus.chromosome.fasta.sa
    7) Example:
    SVprotocol_110_120_130.sh --sample AHRm.725
                            --config example.config
                            --thread 2
    "
}


OPTIONS=hs:c:j:f
LONGOPTIONS=help,sample:,config:,thread:,force
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [ $? != 0 ]; then
    echo "ERROR: Terminating..."
    exit 2
fi
eval set -- "${PARSED}"

FLAG_sample=0
FLAG_config=0
thread=5
force="n"

while true
do
    case "$1" in
        -s|--sample)
            sample=$2
            FLAG_sample=1
            shift 2
            ;;
        -c|--config)
            config=$2
            FLAG_config=1
            shift 2
            ;;
        -j|--thread)
            thread=$2
            shift 2
            ;;
        -f|--force)
            force="y"
            shift
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

#===============================================================================
#****/ Main 
#===============================================================================
source $config
cd $your_analysis_dir

echo "SVprotocol_part1.sh"
printf "\t--sample\t%s\n" "$sample"
printf "\t--config\t%s\n" "$config"
printf "\t--thread\t%s\n" "$thread"
printf "\t--force\t%s\n" "$force"
printf "%-20s\t%-s\n" "Hostname:" "`hostname`"
printf "%-20s\t%-s\n" "Working directory:" "`pwd`"
printf "INFO\t%-s\t%-s\n" "`date`" "Start time"

#****/ Check required files
if [ ! -f ${ref_genome_fasta}.ann -o ! -f ${ref_genome_fasta}.fai ]
then
    echo "Please index your fasta using bwa with command: bwa index ref.fasta"
	exit 1
fi

if [ ! -f ${ref_genome_fasta/.fasta/}.dict -o ! -f ${ref_genome_fasta/.gz/}.dict ]
then
    echo "Please create dictionary for genome fasta: picard.jar CreateSequenceDictionary "
	exit 1
fi

if [ "$force" == "y" ]
then
    rm -rf 110.raw.reads.mapping/$sample >/dev/null
    rm -rf 120.remove.duplicates/$sample >/dev/null
    rm -rf 130.sv.detection/Manta/$sample >/dev/null
    rm -rf 130.sv.detection/Lumpy/$sample >/dev/null
    rm -rf 130.sv.detection/SVseq2/$sample >/dev/null
    rm -rf 130.sv.detection/svimmer_1sample/$sample >/dev/null
    printf "INFO\t%-s\t%-s\n" "`date`" "Forced delete existed files, including 110.raw.reads.mapping/$sample, 120.remove.duplicates/$sample, 130.sv.detection/Manta/$sample, 130.sv.detection/Lumpy/$sample, 130.sv.detection/SVseq2/$sample, 130.sv.detection/svimmer_1sample/$sample"
fi

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

#************/ Create directory
mkdir -p 110.raw.reads.mapping/$sample >/dev/null
mkdir -p 120.remove.duplicates/$sample >/dev/null
mkdir -p 130.sv.detection/Manta/$sample >/dev/null
mkdir -p 130.sv.detection/Lumpy/$sample >/dev/null
mkdir -p 130.sv.detection/SVseq2/$sample >/dev/null
mkdir -p 130.sv.detection/svimmer_1sample/$sample >/dev/null


#************/ 110.raw.reads.mapping /***********************************
existed_file=$your_analysis_dir/110.raw.reads.mapping/$sample/$sample.bam
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    #* 111.mapping (9h30min) / 12 CPUs
    $bwa mem -t $thread -B 4 -O 6 -E 1 -M -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPU:$sample\tPL:ILLUMINA" $ref_genome_fasta fastq/$sample/${sample}_R1.fastq.gz fastq/$sample/${sample}_R2.fastq.gz | $samtools view -@ $thread -bS - > 110.raw.reads.mapping/$sample/$sample.bam
    if [ $? -ne 0 ]
    then
        printf "ERROR%-s\t%-s\n" "`date`" "Failed: bwa mem $sample"
        exit 1
    fi

    #* 112.Sort bam (20min)
    #* -m INT:    Approximately the maximum required memory per thread, specified either in bytes or with a K, M, or G suffix. [768 MiB];To prevent sort from creating a huge number of temporary files, it enforces a minimum value of 1M for this setting.
    $samtools sort -@ $thread -m 1G -o 110.raw.reads.mapping/$sample/$sample.sorted.bam 110.raw.reads.mapping/$sample/$sample.bam && mv 110.raw.reads.mapping/$sample/$sample.sorted.bam 110.raw.reads.mapping/$sample/$sample.bam
    #* 113.Create index of bam file (1min)
    $samtools index -@ $thread 110.raw.reads.mapping/$sample/$sample.bam
    if [ $? != 0 ]
    then
        printf "%-20s\t%-s\n" "ERROR:" "samtools sort & index $sample"
        exit 1
    fi
    printf "INFO\t%-s\t%-s\n" "`date`" "BWA finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "BWA has been executed before, we skip this step."
fi

#************/ 120: Remove duplicates (2h) /***********************************
cd $your_analysis_dir/120.remove.duplicates
existed_file=$your_analysis_dir/120.remove.duplicates/$sample/$sample.dedup.bam
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    $gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$thread" MarkDuplicates \
    -I ../110.raw.reads.mapping/$sample/$sample.bam \
    -O $sample/$sample.dedup.bam \
    --REMOVE_DUPLICATES false \
    -M $sample/metric.$sample.txt \
    --ASSUME_SORT_ORDER coordinate \
    --COMPRESSION_LEVEL 5 \
    --VALIDATION_STRINGENCY SILENT \
    --CREATE_INDEX true \
    --TMP_DIR $sample/tmp/
    #** remove file
    #** $?: the return value of last command，success is 0, failed is 1
    if [ $? == 0 ]
    then
        rm -rf $your_analysis_dir/120.remove.duplicates/$sample/tmp 2>/dev/null
    fi

    printf "INFO\t%-s\t%-s\n" "`date`" "Remove Duplicated finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Remove Duplicated has been executed before, we skip this step."
fi

#************/ 130: SV detection /***********************************
#****/ 131: Manta 
manta_analysis_path=$your_analysis_dir/130.sv.detection/Manta/$sample
cd $manta_analysis_path
existed_file=$manta_analysis_path/results/variants/diploidSV.vcf.gz
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    manta_input_bam=${your_analysis_dir}/120.remove.duplicates/$sample/${sample}.dedup.bam
    #* 131.1: Configuration Manta
    #* Single Diploid Sample Analysis
    ${manta_install_path}/bin/configManta.py \
    --bam $manta_input_bam \
    --referenceFasta $ref_genome_fasta \
    --runDir ${manta_analysis_path}
    # Successfully created workflow run script.
    # To execute the workflow, run the following script and set appropriate options:
    # /picb/humpopg-bigdata5/liuqi/SVprotocol/130.sv.detection/Manta/AHRm.725/runWorkflow.py
    #* 131.2: workflow execution (1h)
    #! Please use python2 to run Manta
    manta_log=$your_analysis_dir/130.sv.detection/Manta/$sample/${sample}.Manta.log
    $python2 ${manta_analysis_path}/runWorkflow.py -j $thread 1>$manta_log 2>&1 &
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Manta has been executed before, we skip this step."
fi

#****/ 132: Lumpy 
#* Lumpy can't use thread
lumpy_analysis_path=$your_analysis_dir/130.sv.detection/Lumpy/$sample
cd $lumpy_analysis_path
existed_file=$lumpy_analysis_path/${sample}.Lumpy.vcf
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    #**** 132.1: Extract the discordant paired-end alignments (20min)
    lumpy_input_bam=${your_analysis_dir}/120.remove.duplicates/$sample/${sample}.dedup.bam
    $samtools view -@ $thread -b -F 1294 $lumpy_input_bam > ${sample}.discordants.unsorted.bam
    #**** 132.2: Extract the split-read alignments (20min)
    $samtools view -h $lumpy_input_bam | $python2 ${lumpy_install_path}/scripts/extractSplitReads_BwaMem -i stdin | $samtools view -@ $thread -Sb - > ${sample}.splitters.unsorted.bam
    #**** 132.3: Sort both alignments (30min)
    $samtools sort -@ $thread ${sample}.discordants.unsorted.bam -o ${sample}.discordants.bam
    $samtools sort -@ $thread ${sample}.splitters.unsorted.bam -o ${sample}.splitters.bam
    #**** 132.4: Run LUMPY Express (44min)
    lumpy_log=$your_analysis_dir/130.sv.detection/Lumpy/$sample/${sample}.Lumpy.log
    $python2 ${lumpy_install_path}/bin/lumpyexpress \
    -B $lumpy_input_bam \
    -S ${sample}.splitters.bam \
    -D ${sample}.discordants.bam \
    -o ${sample}.Lumpy.vcf 1>$lumpy_log 2>&1 &
    #**** Wait for the above program to run
    wait
    printf "INFO\t%-s\t%-s\n" "`date`" "Manta&Lumpy finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Lumpy has been executed before, we skip this step."
fi

#****/ 133: SVseq2 
svseq2_input_bam=${your_analysis_dir}/120.remove.duplicates/$sample/${sample}.dedup.bam
svseq2_analysis_path=$your_analysis_dir/130.sv.detection/SVseq2/$sample
cd $svseq2_analysis_path
existed_file=$svseq2_analysis_path/${sample}.SVseq2.vcf
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    #**** 133.1: For each chromosome (8h)
    #**** We extract chromosome name from genome.fasta.fai
    #**** GWHAMMN00000001	325144201	17	60	61
    genome_fai=${ref_genome_fasta}.fai
    rm -f tmp_multiRun.sh >/dev/null 2>&1
    minMQ=20
    cat $genome_fai | while read line
    do
    arr=($line)
    chr=${arr[0]}
    echo "$samtools view -@ $thread -bS -q $minMQ -o ${sample}.$chr.bam $svseq2_input_bam $chr:300
    $samtools index -@ $thread ${sample}.$chr.bam
    $svseq2 -r $ref_genome_fasta -c $chr -b ${sample}.$chr.bam --c 3 --o ${sample}.$chr.svseq2.DEL.out > ${sample}.$chr.logfile " >${sample}.$chr.sh
    echo "sh ${sample}.$chr.sh" >>tmp_multiRun.sh
    done
    sh ${this_protocol_script_path}/shell_multipleRun.sh tmp_multiRun.sh $thread

    #**** Wait for the above program to run
    wait
    #** Check to SVseq2 out
    cat $genome_fai | while read line
    do
    arr=($line)
    chr=${arr[0]}
    if [ ! -f ${sample}.${chr}.svseq2.DEL.out ]
    then
        printf "ERROR\t%-s\t%-s\n" "`date`" "SVseq2 failed, Please rerun SVseq2 in chr:$chr"
        exit 1
    elif [[ `ls -l ${sample}.${chr}.svseq2.DEL.out |awk '{print $5}'` == 0 ]]
    then
        printf "ERROR\t%-s\t%-s\n" "`date`" "The size of result of SVseq2 is 0, Please rerun SVseq2 in chr:$chr"
        exit 1
    else
        printf "INFO\t%-s\t%-s\n" "`date`" "SVseq2 successed in chr:$chr"
    fi
    done
    #****/ Merging results in all chromosomes and convert the result to a VCF file
    cd $svseq2_analysis_path
    ls *DEL.out > del.out.list
    $perl ${this_protocol_script_path}/code_convertSVseq22vcf.pl -i del.out.list -o ${sample} -s ${sample}
    printf "INFO\t%-s\t%-s\n" "`date`" "SVseq2 finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "SVseq2 has been executed before, we skip this step."
fi

#************/ Sort VCF /***************************************************
#****/ 131: Manta
cd $manta_analysis_path
existed_file=$manta_analysis_path/${sample}.chr.sort.Manta.vcf.gz
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    manta_out_vcf=$manta_analysis_path/results/variants/diploidSV.vcf.gz
    #* 1- Inversion in Manta should be convert format: 
    #* https://github.com/Illumina/manta/tree/master/src/python/libexec/convertInversion.py
    manta_out_convertINV_vcf=$manta_analysis_path/${sample}.Manta.vcf.gz
    $python2 $code_convertInversion $samtools $ref_genome_fasta $manta_out_vcf | $bgzip -c > $manta_out_convertINV_vcf && $tabix -p vcf -f $manta_out_convertINV_vcf 

    #* Sort Manta result
    genome_dict=${ref_genome_fasta/.fasta/}.dict
    $java -jar $picard SortVcf \
    I=$manta_out_convertINV_vcf \
    O=${sample}.chr.sort.Manta.vcf \
    SEQUENCE_DICTIONARY=$genome_dict \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=false
    #**** Compressed VCF and create index
    $bgzip ${sample}.chr.sort.Manta.vcf && $tabix -p vcf -f ${sample}.chr.sort.Manta.vcf.gz

    printf "INFO\t%-s\t%-s\n" "`date`" "Sort Manta finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Sort Manta has been executed before, we skip this step."
fi

#****/ 132: Sort Lumpy
cd $lumpy_analysis_path
existed_file=$lumpy_analysis_path/${sample}.chr.sort.Lumpy.vcf.gz
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    genome_dict=${ref_genome_fasta/.fasta/}.dict
    $java -jar $picard SortVcf \
    I=${sample}.Lumpy.vcf \
    O=${sample}.chr.sort.Lumpy.vcf \
    SEQUENCE_DICTIONARY=$genome_dict \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=false
    #**** Compressed VCF and create index
    $bgzip ${sample}.chr.sort.Lumpy.vcf && $tabix -p vcf -f ${sample}.chr.sort.Lumpy.vcf.gz
    #**** remove temporary files
    rm -f ${sample}.splitters.* ${sample}.discordants.* >/dev/null 2>&1

    printf "INFO\t%-s\t%-s\n" "`date`" "Sort Lumpy finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Sort Lumpy has been executed before, we skip this step."
fi

#****/ 133: Sort SVseq2
cd $svseq2_analysis_path
existed_file=$svseq2_analysis_path/${sample}.chr.sort.SVseq2.vcf.gz
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    genome_dict=${ref_genome_fasta/.fasta/}.dict
    $java -jar $picard SortVcf \
    I=${sample}.SVseq2.vcf \
    O=${sample}.chr.sort.SVseq2.vcf \
    SEQUENCE_DICTIONARY=$genome_dict \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=false
    #**** Compressed VCF and create index
    $bgzip ${sample}.chr.sort.SVseq2.vcf && $tabix -p vcf -f ${sample}.chr.sort.SVseq2.vcf.gz
    #**** remove temporary files
    rm -f tmp_multiRun.sh ${sample}.*.sh >/dev/null 2>&1
    rm ${sample}.*.bam* >/dev/null 2>&1
    rm -f del.out.list *DEL.out >/dev/null 2>&1

    printf "INFO\t%-s\t%-s\n" "`date`" "Sort SVseq2 finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "Sort SVseq2 has been executed before, we skip this step."
fi

#************/ Filtering SVs with sizes larger than 50bp and smaller than 2Mb (1min)
svimmer1_analysis_path=$your_analysis_dir/130.sv.detection/svimmer_1sample/$sample
mkdir -p $svimmer1_analysis_path
cd $svimmer1_analysis_path
rm -f vcf.list >/dev/null 2>&1
for software in Manta Lumpy SVseq2
do
ln -s $your_analysis_dir/130.sv.detection/$software/$sample/${sample}.chr.sort.${software}.vcf.gz* ./ >/dev/null 2>&1
$perl ${this_protocol_script_path}/code_SV_filtering.pl --tool $software --vcf ${sample}.chr.sort.${software}.vcf.gz --len 50 --xlen 2000000
$tabix -p vcf -f ${sample}.chr.sort.${software}.filtered.vcf.gz
echo "${sample}.chr.sort.${software}.filtered.vcf.gz" >> vcf.list
done

#**** Merging SVs discovered by Manta, Lumpy, and SVseq2 using svimmer software
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
$python3_9 $svimmer --threads $thread vcf.list $chr | $bgzip -c > ${sample}.chr.svimmer.vcf.gz

printf "INFO\t%-s\t%-s\n" "`date`" "svimmer 1 sample finished"

#**** Filtering SVs with sizes larger than 50bp and smaller than 2Mb
$perl ${this_protocol_script_path}/code_SV_filtering.pl --tool svimmer --vcf ${sample}.chr.svimmer.vcf.gz --len 50 --xlen 2000000
$tabix -p vcf -f ${sample}.chr.svimmer.filtered.vcf.gz

#**** remove temporary files
# rm -f ${sample}.chr.sort.*.filtered.vcf.gz* ${sample}.chr.svimmer.vcf.gz* vcf.list >/dev/null 2>&1
# for software in Manta Lumpy SVseq2
# do
# rm -f ${sample}.chr.sort.${software}.vcf.gz* >/dev/null 2>&1
# done
printf "INFO\t%-s\t%-s\n" "`date`" "SV filtering for one sample finished"

#************/ Ending
printf "INFO\t%-s\t%-s\n" "`date`" "Congratulations"
printf "INFO\t%-s\t%-s\n" "`date`" "End"

