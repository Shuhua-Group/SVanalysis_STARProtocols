#!/usr/bin/Rscript

rm(list=ls())

options(stringsAsFactors = FALSE,scipen=200)
#===========================================================================
# Package
# Make sure you have installed following packages in your system
#===========================================================================
# load package
suppressMessages(library(dplyr,quietly = TRUE))
suppressMessages(library(edgeR,quietly = TRUE))

command=matrix(c( 
  'help', 'h', 0,'logical', 'Usage',
  'dir', 'd', 0,'character', 'the directory include inputfile and outputfile',
  'htseq_mat', 'm', 1, 'character', 'a gene count mat',
  'cds_bed', 'b', 1, 'character', 'cds bed file',
  'prefix', 'p', 1, 'character', 'prefix of files'),byrow=T,ncol=5)
args=getopt(command)

if (!is.null(args$help) || is.null(args$dir) || is.null(args$htseq_mat) || 
    is.null(args$cds_bed) || is.null(args$prefix))
{
  cat(paste(getopt(command, usage = T), "\n"))
  cat("Note:\n")
  cat("htseq_mat's format:\n")
  cat("  Gene\tSymbol\tsample1_count\tsample2_count...\tsampleN_count\n")
  cat("cds_bed's format:\n")
  cat("  chr\tstart\tend\tgeneName:feature\tfeature\n")
  cat("\nExample:\n")
  cat("Rscript code_convertHTSeq2GFOLD.r --dir /140.gfold.project --htseq_mat project.htseq.gene.count.mat --cds_bed /path/ref_genome_prefix_name.gff.cds.bed --prefix RM\n")
  q(save="no",status=1)
  
}

dir<-args$dir
htseq_mat<-args$htseq_mat
cds_bed<-args$cds_bed
prefix<-args$prefix


lapply_exonSum<-function(gene1,data)
{
    data.cds.tmp<-data %>% filter(gene==gene1) %>% arrange(V2)
        
    min_pos<-min(data.cds.tmp$V2)
    max_pos<-max(data.cds.tmp$V3)
    gene_length<-max_pos-min_pos+1
    
    data.cds.tmp<-data.cds.tmp %>%
        mutate(start_new=V2-min_pos+1,end_new=V3-min_pos+1)
    tmp_matrix<-matrix(0,nrow=gene_length,ncol=1)
    for(i in 1:nrow(data.cds.tmp))
    {
        tmp_matrix[data.cds.tmp$start_new[i]:data.cds.tmp$end_new[i]]<-1
    }
    exon_sum<-colSums(tmp_matrix)
    tmp.summary<-data.frame(gene=gene1,exonSum=exon_sum)
    return(tmp.summary)
}
get_gene_sumOfAllExon <- function(cds.file)
{
    ## 0-based bed
    data.cds<-read.csv(cds.file,sep="\t",header = FALSE)
    ## 1-based bed
    data.cds$V2<-data.cds$V2+1
    data.cds <- data.cds %>%
        mutate(gene=sapply(str_split(V4, ":"), function(x) x[1])) %>%
        unique()
    ## The length sum of all the exons of this gene
    gene_list<-unique(data.cds$gene)
    data.cds.summary<-lapply(gene_list, lapply_exonSum,data=data.cds)
    data.cds.summary2<- do.call(rbind,data.cds.summary)

    return(data.cds.summary2)
}

get_Gfold_input<-function(
        count_file="metaInfect.htseq.gene.count.mat",
        exonSum=exonSum.HL.1,
        prefix=out
    )
{
    ## read count file
    data.count<-read.csv(count_file,sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
    sample_index<-c(3:ncol(data.count))
    ## ====CPM==========================================================
    ## count per million read (normalized count)
    ## Compute counts per million (CPM)
    ## read_count * 10^6 / (lib.size * norm.factors)
    ## The effective library sizes for use in downstream analysis are lib.size * norm.factors where lib.size contains the original library sizes and norm.factors is the vector of scaling factors computed by this function
    ## log = log2
    y <- DGEList(counts=data.count[,3:ncol(data.count)], genes=data.count[,1:2])
    ## filtering gene with lower expression
    y$samples$group<-c(1:nrow(y$samples))
    keep <- filterByExpr(y)
    y <- y[keep, ,keep.lib.sizes=FALSE]
    # 
    y.CPM<-cpm(y,log = FALSE,normalized.lib.sizes = TRUE)
    data.CPM <- data.frame(y$genes,y.CPM)
    
    data.count<-data.count[keep,]
    
    ## ===FPKM=========================================================
    ## + exonSum
    data.CPM <- data.CPM %>%
        left_join(.,exonSum,by=c("Symbol"="gene"))
    exon_index<-ncol(data.CPM)
    for(i in sample_index)
    {
        data.gene.i<-data.count[,1:2]
        data.count.i<-data.count[,i]
        data.CPM.i<- data.CPM[,i]
        data.exonSum<-data.CPM[,exon_index]
        ## FPKM
        data.FPKM<-data.CPM.i * (10^3) / data.exonSum
        ## sample out
        sample.name<-colnames(data.count)[i]
        out.file.name<-paste(sample.name,".gfold.input",sep="")
        out.data<-data.frame(data.gene.i,count=data.count.i,
                             exonLen=data.exonSum,FPKM=data.FPKM)
        write.table(x=out.data,
                  file = out.file.name,
                  col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
    }
    ## ===TMM-CPM========================================================
    ## TMM 
    y <- calcNormFactors(y)
    y.CPM.TMM<-cpm(y,log = FALSE)
    y.CPM.TMM.out <- data.frame(y$genes,y.CPM.TMM)
    out.file.name<-paste(prefix,".CPM_TMMnorm.out",sep="")
    write.table(x=y.CPM.TMM.out,file = out.file.name,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
}

exonSum.file<-paste(prefix,".exonSum.RData",sep = "")
if(file.exists(exonSum.file))
{
    load(exonSum.file)
}else
{
    exonSum.RM<-get_gene_sumOfAllExon(cds.file = cds_bed)
    save(exonSum.RM,file = exonSum.file)
}
get_Gfold_input(
        count_file=htseq_mat,
        exonSum=exonSum.RM,
        prefix=prefix
    )

