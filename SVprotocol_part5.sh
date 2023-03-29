#!/usr/bin/sh

#===============================================================================
#* This pipeline is used to conduct differential gene expression analysis
#===============================================================================

function usage()
{
    echo "
******This pipeline is used to conduct differential gene expression (DGE) analysis
Usage: SVprotocol_part5.sh   <The long option is recommended>
    -s, --samplelist  (required) sample list, one sample per line
    -c, --config      (required) a config file including basic setting.
    -p, --project     (required) a project name as the suffix of subdirectory. default = project
    -j, --thread      (optional) the number of thread, default = 5
    -h, --help        (optional) This small usage guide

    1) In your_gene_expression_analysis_dir, you must have a subdirectory named as fastq/<sample_name> including two FASTQ file (<sample_name>_1.fastq.gz and <sample_name>_2.fastq.gz)
    2) Please conduct quality control by yourself. 
        This pipeline directly use the fastq files with no quality control.
    3) Please prepare a config file following the example file
        1. Containing the path of software to be used
        2. You must use the variable name shown in the example file because, in this pipeline, we will use these variable names
    4) Please just include chromosome-level sequence or sequence with size 
        larger than 250kb in the reference genome fasta file
    5) <project>: Usually, we will provide a series of samples and run DGE for 
        one destination, e.g. normal vs infected samples. For the convenience of subsequent analysis, this pipeline will put all results into a subdirectory with the suffix name <project>.
    6) Example:
    SVprotocol_part5.sh --samplelist SRR.list
                            --config example.config
                            --project RM
                            --thread 5
    "
}


OPTIONS=hs:c:j:
LONGOPTIONS=help,samplelist:,config:,thread:
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [ $? != 0 ]; then
    echo "ERROR: Terminating..."
    exit 2
fi
eval set -- "${PARSED}"

FLAG_sample=0
FLAG_config=0
thread=5
project="project"

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
        -j|--thread)
            thread=$2
            shift 2
            ;;
        -p|--project)
            project=$2
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
config=/picb/humpopg-bigdata5/liuqi/SVprotocol/example2.config
sample=AHRm.725
thread=12

source $config
cd $your_gene_expression_analysis_dir

echo "SVprotocol_part5.sh"
printf "\t--samplelist\t%s\n" "$samplelist"
printf "\t--config\t%s\n" "$config"
printf "\t--project\t%s\n" "$project"
printf "\t--thread\t%s\n" "$thread"
printf "%-20s\t%-s\n" "Hostname:" "`hostname`"
printf "%-20s\t%-s\n" "Working directory:" "`pwd`"
printf "INFO\t%-s\t%-s\n" "`date`" "Start time"

#************/ Create directory
mkdir -p $your_gene_expression_analysis_dir/110.hisat2.index_genome 2>/dev/null
mkdir -p $your_gene_expression_analysis_dir/120.hisat2.genome_out.${project} 2>/dev/null
mkdir -p $your_gene_expression_analysis_dir/130.htseq.genome_out.${project} 2>/dev/null
mkdir -p $your_gene_expression_analysis_dir/140.gfold.${project} 2>/dev/null


#************/ 110.hisat2.index_genome /***************************
hisat2_index_dir=$your_gene_expression_analysis_dir/110.hisat2.index_genome
hisat2_index_prefix=$ref_genome_prefix_name

cd $hisat2_index_dir

existed_file=$hisat2_index_dir/${hisat2_index_prefix}.1.ht2
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    #* build index only with genome file (1h) / 3 CPUs
    $hisat2_build -p $thread $ref_genome_fasta $hisat2_index_dir/$hisat2_index_prefix
    printf "INFO\t%-s\t%-s\n" "`date`" "hisat2 index with genome finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "hisat2 index has been executed before, we skip this step."
fi

#************/ 120.hisat2.genome_out /***********************************
cd $your_gene_expression_analysis_dir
hisat2_analysis_dir=$your_gene_expression_analysis_dir/120.hisat2.genome_out.${project}
genome_index=$your_gene_expression_analysis_dir/110.hisat2.index_genome/$hisat2_index_prefix

cat $samplelist | while read sample
do
    INDIR=$your_gene_expression_analysis_dir/fastq/$sample
    OUTDIR=$hisat2_analysis_dir
    existed_file=$hisat2_analysis_dir/$sample.sorted.bam
    check_existed_file $existed_file
    if [ $? -eq 1 ]
    then
        ## run hisat2
        $hisat2 -p $thread -x $genome_index -1 $INDIR/${sample}_1.fastq.gz -2 $INDIR/${sample}_2.fastq.gz --seed 110 -S $OUTDIR/$sample.sam && $samtools view -@ $thread -bS -h $OUTDIR/$sample.sam > $OUTDIR/$sample.bam && $samtools sort -@ $thread -T $sample $OUTDIR/$sample.bam > $OUTDIR/$sample.sorted.bam && $samtools index $OUTDIR/$sample.sorted.bam && rm -f $OUTDIR/$sample.sam $OUTDIR/$sample.bam >/dev/null 2>&1
        printf "INFO\t%-s\t%-s\n" "`date`" "$sample hisat2 mapping finished"
    else
        printf "INFO\t%-s\t%-s\n" "`date`" "$sample hisat2 mapping has been executed before, we skip this step."
    fi
done


#************/ 130.htseq.genome_out /***********************************
cd $your_gene_expression_analysis_dir
htseq_analysis_dir=$your_gene_expression_analysis_dir/130.htseq.genome_out.${project}

cat $samplelist | while read sample
do
    INDIR=$hisat2_analysis_dir
    OUTDIR=$htseq_analysis_dir
    existed_file=$OUTDIR/$sample.sorted.htseq.gene.count
    check_existed_file $existed_file
    if [ $? -eq 1 ]
    then
        ## run HTSeq basing on gene_id
        # -s no  : indicates we're using an unstranded RNA-seq library.
        # -r pos : tells htseq-count that our BAM file is coordinate sorted.
        # -f bam : indicates that our input file is in BAM format.
        # --type FEATURE_TYPE:
        #        Feature type (3rd column in GTF file) to be used, all
        #        features of other type are ignored (default, suitable
        #        for Ensembl GTF files: exon)
        #--idattr IDATTR:
        #       read count是基于feature来统计的;
        #       GTF attribute to be used as feature ID (default,
        #       suitable for Ensembl GTF files: gene_id). All feature
        #       of the right type (see -t option) within the same GTF
        #       attribute will be added together. The typical way of
        #       using this option is to count all exonic reads from
        #       each gene and add the exons but other uses are
        #       possible as well.
        #--nonunique=<nonunique mode> Mode to handle reads that align to or are assigned to more than one feature in the overlap <mode> of choice (see -m option). <nonunique mode> are none, all, fraction, and random (default: none)

        $htseq_count -s no -r pos -f bam --nprocesses=$thread --mode=union --type exon --idattr gene_id --additional-attr=gene_name --nonunique=all $INDIR/$sample.sorted.bam $ref_genome_gtf > $OUTDIR/$sample.sorted.htseq.gene.count

        printf "INFO\t%-s\t%-s\n" "`date`" "$sample htseq count finished"
    else
        printf "INFO\t%-s\t%-s\n" "`date`" "$sample htseq count has been executed before, we skip this step."
    fi
done

#** Combine all sample's gene count
cd $htseq_analysis_dir
mat_input_suffix="sorted.htseq.gene.count"
$perl $this_protocol_script_path/code_convertCount2Matrix.pl --dir $htseq_analysis_dir --srr_list $your_gene_expression_analysis_dir/$samplelist --out_prefix $project --input_suffix $mat_input_suffix
#** output file: <project>.htseq.gene.count.mat


#************/ 140.gfold /***********************************
cd $your_gene_expression_analysis_dir
gfold_analysis_dir=$your_gene_expression_analysis_dir/140.gfold.${project}
cd $gfold_analysis_dir
ln -s $htseq_analysis_dir/$project.htseq.gene.count.mat >/dev/null 2>&1
#**** Convert gff to bed file to be used in analysis
# output files
# out_prefix.gff.cds.bed
# out_prefix.gff.mRNA.bed
# out_prefix.gff.UTR.bed
# out_prefix.gff.promoter.bed
# out_prefix.gff.upstream.bed
# out_prefix.gff.downstream.bed
out_gff_convert_dir=`dirname $ref_genome_gff`
out_prefix=$out_gff_convert_dir/$ref_genome_prefix_name
existed_file=$out_prefix.gff.cds.bed
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    $perl $this_protocol_script_path/code_gff2SVAnnBed.pl $ref_genome_gff $out_prefix
    printf "INFO\t%-s\t%-s\n" "`date`" "gff.cds.bed finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "gff.cds.bed has been executed before, we skip this step."
fi

#**** Convert gene count to gfold input
cd $gfold_analysis_dir

gff_cds_bed=$out_gff_convert_dir/$ref_genome_prefix_name.gff.cds.bed
sample=`head -n1 $your_gene_expression_analysis_dir/$samplelist`
existed_file=$OUTDIR/$sample.gfold.input
check_existed_file $existed_file
if [ $? -eq 1 ]
then
    $Rscript $this_protocol_script_path/code_convertHTSeq2GFOLD.r --dir $your_gene_expression_analysis_dir/140.gfold.$project --htseq_mat $project.htseq.gene.count.mat --cds_bed $gff_cds_bed --prefix $project

    printf "INFO\t%-s\t%-s\n" "`date`" "$sample gfold input data finished"
else
    printf "INFO\t%-s\t%-s\n" "`date`" "$sample gfold input data has been executed before, we skip this step."
fi

#************/ Ending
printf "INFO\t%-s\t%-s\n" "`date`" "Congratulations"
printf "INFO\t%-s\t%-s\n" "`date`" "End"
