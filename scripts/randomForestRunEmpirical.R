library(optparse)

parser <- OptionParser()

parser <- add_option(
    parser,
    "--modelFile",
    type = "character",
    action = "store",
    help = ".RData file with trained model."
)

parser <- add_option(
    parser,
    "--fisFile",
    type = "character",
    action = "store",
    help = "output from plink --het for empirical dataset."
)

parser <- add_option(
    parser,
    "--tajimaFile",
    type = "character",
    action = "store",
    help = "output from vcftools --TajimaD for empirical dataset."
)

parser <- add_option(
    parser,
    "--rohFile",
    type = "character",
    action = "store",
    help = "ROH output file from BCFtools."
)

parser <- add_option(
    parser,
    "--chromosomeLengths",
    type = "character",
    action = "store",
    help = "bp length for each chromosome"
)

parser <- add_option(
    parser,
    "--cpus",
    type = "integer",
    action = "store",
    default = 1,
    help = "Number of CPUs [default %default]"
)

parser <- add_option(
    parser,
    "--out",
    type = "character",
    action = "store",
    default = "selfRates.txt",
    help = "desired output file"
)


args <- parse_args(parser)

#### for testing
## args$modelFile <- "caretmodels/modelPreProc20230327171747JeCod.RData"
## args$modelFile <- "caretmodels/modelPreProc20230327171706uZ2bC.RData"
## args$fisFile <- "empirical/fb_run7_all_mj_qu_id_snpeff_uniq.St.NEU.uniq.AA.vcf.gz.het"
## args$tajimaFile <- "empirical/fb_run7_all_mj_qu_id_snpeff_uniq.St.NEU.uniq.AA.vcf.gz.Tajima.D"
## args$rohFile <- "empirical/fb_run7_all_nn_bi_mj_sn_qu_mm_id_PL.vcf_roh.txt"
## args$chromosomeLengths <- "empirical/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.lengths"
## args$out <- "empirical/selfrates/St.txt"
## Rscript randomForestRunEmpirical.R --modelFile caretmodels/modelPreProc20230327171747JeCod.RData --fisFile empirical/fb_run7_all_mj_qu_id_snpeff_uniq.St.NEU.uniq.AA.vcf.gz.het --tajimaFile empirical/fb_run7_all_mj_qu_id_snpeff_uniq.St.NEU.uniq.AA.vcf.gz.Tajima.D --rohFile empirical/fb_run7_all_nn_bi_mj_sn_qu_mm_id_PL.vcf_roh.txt --chromosomeLengths empirical/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.lengths --out empirical/selfrates/St1.txt

#### -- RUN --
## source('~/code/r/source_me.R')

suppressPackageStartupMessages({
    require(dplyr)
    require(data.table)
    require(caret)
})

## source("../rffunctions.R")
## source("../confMat.R")

setDTthreads(threads = args$cpus)

cat("\n\n==========================================\n\nParameters:\n\n")
cat(paste0(names(args),": ",args,collapse = "\n"),"\n")
cat("\nStarting data prep.\n")

## prep stats
if (!is.null(args$fisFile)) {
    het0 <- fread(args$fisFile,data.table=F,sep=" ")
    het <- het0 %>%
        select(sample=FID,Fis=F) %>%
        mutate(Fis=as.numeric(Fis)) %>% 
        filter(!is.na(Fis)) %>%
        summarise(Fis=mean(Fis))
    idMask <- het0 %>%
        select(sample=FID)
    nInd <- idMask %>% pull %>% unique %>% length
} else {
    het <- NULL
    cat("Did not supply fisFile, trying without it.\n")
    idMask <- NULL
}

if (!is.null(args$tajimaFile)) {
    taj <- fread(args$tajimaFile,data.table=F,sep="\t") %>%
        select(bin=BIN_START,chr=CHROM,TajimaD) %>%
        mutate(TajimaD=as.numeric(TajimaD)) %>%
        filter(!is.na(TajimaD)) %>%
        group_by(chr) %>% 
        summarise(TajimaD=mean(TajimaD))
} else {
    taj <- NULL
    cat("Did not supply tajimaFile, trying without it.\n")
}

stats <- data.frame(taj,het)

cat("Prepared stats.\n")


#### prep ROH data
gsize <- fread(args$chromosomeLengths,data.table=F) %>%
    rename(chr=1,chrlength=2)

roh <- fread(args$rohFile,data.table = F)
names(roh) <- c('rg','sample','chr','start','end','length','snps','quality')

if (!is.null(idMask)) { # make sure only individuals are used that occur in Fis+ROH
    roh <- roh %>%
        filter(sample%in%idMask$sample)
}

roh <- roh %>%
    inner_join(gsize,by='chr') %>%             # get lengths of chromosomes
    right_join(idMask,by='sample') %>%         # merge with all inds, no ROH for indv?
    group_by(sample,chr) %>%                   # individual-wise
    mutate(                                    #
        gap=coalesce(lead(start)-end,0),       # calculate gaps
        proportionInROH=coalesce(sum(length)/chrlength,0), # proportion relative to chromosome length
        roh_count_ind=coalesce(length(sample)/chrlength,0),# counts of ROH
        length=coalesce(length,0)                          # raw length
    ) %>%
    group_by(chr) %>% # population-wise
    summarise(
        lengthVar=coalesce(var(length/chrlength),0),
        countVar=coalesce(var(roh_count_ind),0),
        propVar=coalesce(var(proportionInROH),0),
        roh_count_ind=mean(roh_count_ind),
        proportionInROH=mean(proportionInROH),
        lengthMedian=median(length/chrlength),
        gapMedian=coalesce(median(gap/chrlength,na.rm=T),0),
        gapVar=coalesce(var(gap/chrlength,na.rm=T),0)
    )

cat("Prepared ROHs.\n")

## merge
if (nrow(stats)==0) {
    fullData <- roh
} else {
    fullData <- inner_join(roh,stats,by='chr')
}
dTest <- data.frame(fullData)

## -- PREPARE MODELING --
cat("Loading model.\n")
load(args$modelFile)

modelrf <- model[[1]]
modelInfoStr <- model[[2]]

## -- PREPROCESSING --
cat("Starting preprocessing.\n")
preProcValues <- model[[3]]

dTestTransformed <- predict(preProcValues, dTest)   # new test set
names(dTest) <- paste0(names(dTest), ".orig")
dTest <- bind_cols(dTest,
                   dTestTransformed)

## -- PREDICTION --
## predict with test set and evaluate model
cat("Starting model.\n")
pred <- predict(modelrf,dTest)

## save output
cat("Saving output: ",args$out,"\n")
out <- bind_cols(fullData,selfRate=pred)

fwrite(out,paste0(args$out,".params"))
fwrite(data.frame(pred),args$out,col.names=F)

cat("Done!\n===================================\n===================================\n\n")
