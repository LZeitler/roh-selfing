source('source_me_light.R')


## training/test data - REGRESSION
tstamp <- "20231002191703IfZEj"

load(paste0("output/predict/",tstamp,".dTest.RData"))
predicted <- fread(paste0("output/predict/",tstamp,".predicted.csv"),data.table=F)

dPred <- cbind(dTest,pred=predicted$prediction) %>%
    mutate(absErr=abs(pred-selfRate))



## CORRPLOT
vars <- c(
    "lengthVar", "countVar", "propVar", "roh_count_ind","proportionInROH",
    "lengthMedian", "gapMedian", "gapVar", "TajimaD",  "Fis", "selfRate",
    "absErr")
varnames1 <- c(
    "Variance Lenths"="lengthVar",
    "Variance Counts"="countVar",
    "$ Variance~F[ROH]"="propVar",
    "Counts of ROH"="roh_count_ind",
    "$F[ROH]"="proportionInROH",
    "Median Length of ROH"="lengthMedian",
    "Median non-ROH Gap length"="gapMedian",
    "Variance Gaps"="gapVar",
    "Tajima's D"="TajimaD",
    "Inbreeding coefficient"="Fis",
    "Selfing Rate"="selfRate",
    "Absolute Error"="absErr"
)
dPredC <- dPred %>% rename(all_of(varnames1))
pred.cor <- cor(dPredC[,names(varnames1)])

pdf(paste0("plots/predVsErr.corplot2.",tstamp,".pdf"),width=9,height=9)
corrplot::corrplot(pred.cor,order='AOE',addCoef.col='black',tl.col="black",method='ellipse')
dev.off()

## varnames for ggplot
varnames <- c(
    "Variance Lenths"="lengthVar",
    "Variance Counts"="countVar",
    "Variance~F[ROH]"="propVar",
    "Counts of ROH"="roh_count_ind",
    "F[ROH]"="proportionInROH",
    "Median Length of ROH"="lengthMedian",
    "Median non-ROH Gap length"="gapMedian",
    "Variance Gaps"="gapVar",
    "Tajima's D"="TajimaD",
    "Inbreeding coefficient"="Fis",
    "Selfing Rate"="selfRate",
    "Absolute Error"="absErr"
)



## FIGURE 1; MATRIX PROP VS COUNT VS SELFING VS DEMOGRAPHY
transformIf <- function(data,variable) {
    if (any(data[,variable]<0)
        ) {
        mytransformx <- mytransformy <- NULL
    } else {
        mytransformx <- scale_x_log10()
        mytransformy <- scale_y_log10()
    }
    return(list(mytransformx,mytransformy))
}

plotMatrix <- function(data,ix,iy,shape="N_extra") {
    data <- data %>% 
        mutate(e=ifelse(!demo_3class%in%c("bneck","onepop"),demo,demo_3class))
    adlvl <- data %>% pull(e) %>% unique
    lvls <- c("Constant","Bottleneck","Admixture","Isolation",adlvl)
    cat(lvls)
              
    data <- data %>%
        mutate(
            demo_extra=
                factor(case_match(e,
                           "onepop"~"Constant",
                           "bneck"~"Bottleneck",
                           "twoPop"~"Admixture",
                           "twoPopIsol"~"Isolation",
                           .default = e),
                levels=lvls),
            selfRate_extra=paste("Selfing Rate:",selfRate*100,"%"),
            roh_count_ind.absolute.orig=roh_count_ind.orig*1e7,
            N_extra=paste0("N=",N))

    mytransformx <- transformIf(data,ix)[[1]] # transform variable if data allows for it
    mytransformy <- transformIf(data,iy)[[2]] # same for y axis
    
    p4 <- ggplot(data)+
        ## scale_color_viridis_c(
            ## option = "plasma",
            ## direction = -1,
            ## end = .95)+
        background_grid(minor="x")+
        geom_point(
            aes(
                !!sym(ix),
                !!sym(iy),
                shape=!!sym(shape)
                ## color=selfRate
            ),
            alpha=.75,
            size=3
        )+
        scale_shape_manual(values=c(3,4,1,2))+
        ## mytransformx+
        mytransformy+
        ggh4x::facet_grid2(selfRate_extra~demo_extra,
                           axes="all",
                           remove_labels="all")
    return(p4)
}

p0Data <- filter(
    dPred
   ,(N==1000&demo_3class=='bneck'&subpop=='contr') |  # any contracted bottleneck small
    (N==10&demo=='onepop') |                          # onepop tiny
    (N==100&demo=='onepop') |                         # onepop smaller
    (N==1000&demo=='onepop') |                        # onepop small
    (N==10000&demo=='onepop') |                       # onepop big
    (N==1000&demo=="twoPop"&subpop=="EarlyAdmixP1") | # admixture small early
    (N==1000&demo=="twoPopIsol")                      # any isolation small
   ,selfRate%in%c(0, 0.5, .99))

p0Data$demo_extra

p0 <- plotMatrix(p0Data,
    "proportionInROH.orig","roh_count_ind.absolute.orig") +
    labs(x=parse(text=names(varnames[varnames=="proportionInROH"])),
         y=names(varnames[varnames=="roh_count_ind"]))+
    theme(legend.position=c(.09,.9),
          legend.title = element_blank(),
          legend.background = element_rect(fill=alpha("grey70",.5)),
          legend.margin = margin(t = -4, r = 3, b = 2, l = 3, unit = "pt"))+
    scale_x_continuous(breaks=c(0,.5,1),labels=c(0,0.5,1))
p0

## demographic comic - load from inkscape
co_admix <- magick::image_read("comics/admix1.png")
co_bot <- magick::image_read("comics/bottleneck1.png")
co_const <- magick::image_read("comics/constant1.png")
co_isol <- magick::image_read("comics/isolation1.png")
co_arrow <- magick::image_read("comics/arrow.png")

p_admix <- ggdraw()+draw_image(co_admix)
p_bot <- ggdraw()+draw_image(co_bot)
p_const <- ggdraw()+draw_image(co_const)
p_isol <- ggdraw()+draw_image(co_isol)
p_arrow <- ggdraw()+draw_image(co_arrow)
                    

## assemble with patchwork
require(patchwork)

f1 <- (p_const|p_bot|p_admix|p_isol)/
    p0+
    plot_layout(heights = c(.8, 2))
f1 <- f1+
    inset_element(
        p_arrow,
        left = 0, bottom = 1, right = .1, top = 1.35,
        align_to = 'full'
    )
f1

saveplot(f1,"plots/20240220_f1",8,8)


## SUPPLEMENT: EXTRA BOTTLENECK TYPES PROP VS LENGTH
p0DataExtra <- filter(
    dPred
   ,(N==1000&demo_3class=='bneck')
   ,selfRate%in%c(0, 0.5, .99)) %>%
    mutate(demo.orig=ifelse(demo.orig=="bneckGrad","gradual",
                     ifelse(demo.orig=="bneckLong","long",
                     ifelse(demo.orig=="bneckShort","short",NA))),
           demo_3class=ifelse(subpop=="contr","Contraction","Recovery"),
           demo=ifelse(subpop=="contr","Contraction","Recovery"))

p0 <- plotMatrix(
    data=
        p0DataExtra,
    ix=
        "proportionInROH.orig",
    iy=
        "roh_count_ind.absolute.orig",
    shape=
        "demo.orig")+
    labs(x=parse(text=eval(names(varnames[varnames=="proportionInROH"]))),
         y=names(varnames[varnames=="roh_count_ind"]),
         shape="Bottleneck\nTypes")
p0

saveplot(p0,"plots/20240220_f1.png_plotMatrix_BneckTypesBeforeAfter",8,8)

## N=10000
p0DataExtra10k <- filter(
    dPred
   ,(N==10000&demo_3class=='bneck')
   ,selfRate%in%c(0, 0.5, .99)) %>%
    mutate(demo.orig=ifelse(demo.orig=="bneckGrad","gradual",
                     ifelse(demo.orig=="bneckLong","long",
                     ifelse(demo.orig=="bneckShort","short",NA))),
           demo_3class=ifelse(subpop=="contr","Contraction","Recovery"),
           demo=ifelse(subpop=="contr","Contraction","Recovery"))

p0 <- plotMatrix(
    data=
        p0DataExtra10k,
    ix=
        "proportionInROH.orig",
    iy=
        "roh_count_ind.absolute.orig",
    shape=
        "demo.orig")+
    labs(x=parse(text=eval(names(varnames[varnames=="proportionInROH"]))),
         y=names(varnames[varnames=="roh_count_ind"]),
         shape="Bottleneck\nTypes")
p0

saveplot(p0,"plots/20240220_f1.png_plotMatrix_BneckTypesBeforeAfterN10000",8,8)


## FIGURE 2; COMIC, PREDICTIONS, VARIMP
## functions
p_viol_plotter <- function(data) {
    p <- ggplot(data)+
        geom_violin(aes(x=paste(data),y=prediction))+
        geom_hline(aes(yintercept=data),color='red')+
        guides(x=guide_axis(angle=45))+
        labs(x="True Selfing",y="Prediction")+
        facet_grid(.~paste(data),scales="free",space="free")+
        theme(strip.text = element_blank(),
              panel.spacing = unit(0,"lines"))+
        background_grid(major="y")
    return(p)
}
conf_matrix <- function(cm){
    cm_d <- as.data.frame(cm$table)
    cm_st <-data.frame(cm$overall)
    cm_st$cm.overall <- round(cm_st$cm.overall,2)
    cm_p <- as.data.frame(prop.table(cm$table))
    cm_d$Perc <- round(cm_p$Freq*100,2)
    cm_d <- cm_d %>%
        group_by(Reference) %>%
        mutate(rsum=sum(Freq), 
               GroupRatio=Freq/rsum)
    p <-  ggplot(data = cm_d, aes(x = Reference , y =  Prediction, fill = GroupRatio))+
        geom_tile() +
        geom_text(aes(label = Freq), color = 'red', size = 5) +
        scale_fill_distiller(palette="Blues",direction = 1)+
        guides(fill="none")+
        guides(x=guide_axis(angle=45))+
        labs(x="True Demography",y="Prediction")
    return(p)
}
demo3Function <- function(data){
    ifelse(data=="bneck","Bottleneck",
    ifelse(data=="onepop","Constant",
    ifelse(data=="twoPop","Admix./Isol.",NA)))
}
convertDemoStrings <- function(data) {
    strs0 <- data.frame(do.call(rbind,stringr::str_split(data,'_')))
    strs1 <- data.frame(demo=do.call(rbind,stringr::str_split(strs0[,1],'\\d+'))[,1],
                        N=do.call(rbind,stringr::str_split(strs0[,1],'\\D+'))[,2],
                        subpop=strs0[,2])
    strs1$demo_bneckVsOther <- gsub(".*bneck.*","bneck",strs1$demo) # bottnecks lumped
    strs1$demo_3class <- ifelse(grepl("bneck",strs1$demo),"bneck",
                         ifelse(grepl("twoPop",strs1$demo),"twoPop",
                         ifelse(grepl("onepop",strs1$demo),"onepop",NA)))
    strs1$demo_2class <- ifelse(grepl("bneck",strs1$demo),"bneck","const")
    strs1$demo_bneckVsOtherSN <- paste(strs1$demo_bneckVsOther,strs1$subpop,strs1$N)
    strs1$demo_3classSN <- paste(strs1$demo_3class,strs1$subpop,strs1$N)
    strs1$demo_2classSN <- paste(strs1$demo_2class,strs1$subpop,strs1$N)
    strs1$Ncurrent <- paste(ifelse(grepl("contr",strs1$subpop), # current N 
                                   as.numeric(strs1$N)/100,as.numeric(strs1$N)))
    strs1$demo_3class_nice <- ifelse(strs1$demo_3class=="bneck","Bottleneck",
                              ifelse(strs1$demo_3class=="onepop","Constant",
                              ifelse(strs1$demo_3class=="twoPop","Admixture/Isolation",NA)))
    strs1$demo_3class_Ncurrent <- paste(strs1$demo_3class_nice,strs1$Ncurrent)
    strs1$demo_3class_N <- paste(strs1$demo_3class_nice,strs1$N)
    strs1$subpop_nice <- ifelse(grepl("contr",strs1$subpop),"contracted",
                         ifelse(grepl("Early",strs1$subpop),"early","stable"))
    strs1$demo_3class_NS <- paste(strs1$demo_3class_nice,strs1$subpop_nice,strs1$N)
    strs1$demo_3class_NcSi <- interaction(
        factor(strs1$demo_3class_nice,
               levels=c("Constant","Bottleneck","Admixture/Isolation")),
        factor(strs1$subpop_nice,
               levels=c("stable","early","contracted")),
        as.factor(strs1$Ncurrent),drop=T,sep=" ")

    return(strs1)
}

varnamesShort <- c(
    "lengthMedian"="Length",
    "lengthVar"="`Length Var.`",
    "roh_count_ind"="Count",
    "countVar"="`Count Var.`",
    "proportionInROH"="F[ROH]",
    "propVar"="F[ROH]~Var.",
    "gapMedian"="Gap",
    "gapVar"="`Gap Var.`",
    "Fis"="F",
    "TajimaD"="`Tajima's D`",
    "selfRate"="`Selfing Rate`",
    "absErr"="`Abs. Error`"
)

## Prepare Variable Importance
f <- fread("../models/infostrs/20231002_original.txt",data.table = F,header=F)
t <- stringr::str_split(f$V11,",")
vars <- do.call(
    bind_rows,
    sapply(t, function(x) {
        setNames(as.numeric(sapply(stringr::str_split(x,":"),"[[",2)),
                 sapply(stringr::str_split(x,":"),"[[",1))
    }))
var <- bind_cols(vars,f)
var <- pivot_longer(var,names(vars),names_to = "predictor",values_to = "imp")
var$predCount <- paste(stringi::stri_count(var$V4,regex="\\+")+1)
var$predictor <- factor(var$predictor,levels=names(vars))
var$alternate <- var$predictor%in%names(vars)[seq(1,length(names(vars)),2)]
var$model <- factor(var$V4,
                    levels=c(
                        "Fis+TajimaD",
                        "lengthMedian+lengthVar+roh_count_ind+countVar+proportionInROH+propVar+gapMedian+gapVar",
                        "lengthMedian+lengthVar+roh_count_ind+countVar+proportionInROH+propVar+gapMedian+gapVar+Fis+TajimaD"
                    ),
                    labels=c(
                        "RF summary",
                        "RF ROH",
                        "RF full"
                    ))
var <- var %>% filter(!grepl("h1bht",V9)) # ran this twice but only need one

var$perf <- unlist(stringr::str_extract_all(
                                var$V6
                               ,"Accuracy:(\\d+\\.\\d+)|RMSE:(\\d+\\.\\d+)"))
## var$perf_all <- paste0("list(",gsub(":","==",gsub("Rsquared","R^{2}",var$V6)),")")
t <- stringr::str_split_fixed(gsub(":","==",gsub("Rsquared","R^{2}",var$V6)),", ",3)
var$perf1 <- t[,1]
var$perf2 <- t[,2]
var$perf3 <- t[,3]
ncurrentID <- filter(var,V3=="Ncurrent")$V9%>%unique # this is the ID for the sequential model
var <- filter(var, !(V9==ncurrentID & V5=="Classification")) # this is the prediction of Ncurrent

varplotFun <- function(data) {
    data <- data %>% mutate(predictor=factor(predictor,
                                     levels=names(varnamesShort),
                                     labels=varnamesShort))
    p <-     ggplot(data)+
    geom_col(
        aes(predictor,Inf,alpha=alternate),
        position = position_dodge2(.5),
        width = 1,
        show.legend=F
    ) +
    geom_col(
        aes(predictor,imp),
        fill="grey50",
        position=position_dodge2(.5),
        na.rm=T)+
    geom_point(
        aes(predictor,imp,group=model),
        position=position_dodge2(.5),
        na.rm=T,show.legend=F)+
    geom_text(
        aes(10.4,90,label=perf1),
        hjust=1,parse=T,check_overlap=T
    )+
    geom_text(
        aes(10.4,78,label=perf2),
        hjust=1,parse=T,check_overlap=T
    )+
    geom_text(
        aes(10.4,64,label=perf3),
        hjust=1,parse=T,check_overlap=T
    )+
    scale_alpha_manual(values=c(0,.15))+
    guides(x=guide_axis(angle=45))+
    scale_x_discrete(labels=function (x) parse(text = x))+
    background_grid(major="y")+
    labs(x="Predictor",y="Relative Importance",fill="Model")+
    theme(legend.position = "none")
    return(p)
}

## A: comic - outline of workflow
co_workflow <- magick::image_read_pdf("comics/process_comic_right.pdf")
p_workflow <- ggdraw()+draw_image(co_workflow)


## 2B: RF ROH-ONLY TRAIN-TEST, VARIABLE IMPORTANCE
tstamp <- "20231002191703bwG5E"

varT <- var %>%
    filter(
        V9==tstamp,
        !is.na(imp)                                                                    
    )

f2b <- varplotFun(varT)
f2b

## 2C: RF ROH-ONLY TEST-UNSEEN, PREDICITIONS
tstamp <- "20231002191703bwG5E"
tstampNew <- "202310031627367d95b179"

predicted <- fread(
    paste0("output/predict/",tstamp,".run",tstampNew,".predicted.csv"),
    data.table=F)

f2c <- p_viol_plotter(data=predicted)
f2c

## 2D: RF FULL TRAIN-TEST, VARIABLE IMPORTANCE
tstamp <- "20231002191703IfZEj"

varT <- var %>%
    filter(
        V9==tstamp,
        !is.na(imp)                                                                    
    )

f2d <- varplotFun(varT)
f2d

## 2E: RF FULL TEST-UNSEEN, PREDICITIONS
tstamp <- "20231002191703IfZEj"
tstampNew <- "202310031627384c9524cd"

predicted <- fread(
    paste0("output/predict/",tstamp,".run",tstampNew,".predicted.csv"),
    data.table=F)

f2e <- p_viol_plotter(data=predicted)
f2e

## SUPPLEMENT: RF SEQUENTIAL, VARIABLE IMPORTANCE
tstamp <- "202310021917048AtJy"
varT <- var %>%
    mutate(predictor=ifelse(grepl("Ncurrent",as.character(predictor),ignore.case=T),
                            "N",as.character(predictor))) %>%
    group_by(predictor,model,perf1,perf2,perf3,V9) %>%
    summarise(imp=sum(imp)) %>% 
    filter(
        V9==tstamp,
        !is.na(imp)                                                                    
    )
varT$alternate <- varT$predictor%in%names(vars)[seq(1,length(names(vars)),2)]
varT$predictor <- factor(varT$predictor,levels=c("N",names(vars)))


f2s3 <- varplotFun(varT)+scale_x_discrete(labels=parse(text=c("N",varnamesShort)))
f2s3

varT %>% select(predictor,imp) %>% arrange(-imp) # numbers

## SUPPLEMENT: RF SEQUENTIAL, N-CURRENT PREIDICTION (MATRIX)
tstamp <- "N202310021917048AtJy"

predicted <- fread(paste0("output/predict/",tstamp,".predicted.csv"),data.table=F)

predicted_n <- predicted %>%
    mutate_at(c("prediction","data"),as.factor)

cm <- caret::confusionMatrix(predicted_n$prediction,predicted_n$data)

f2s4 <- conf_matrix(cm=cm)
f2s4 <- f2s4+labs(x="True N")
f2s4

## SUPPLEMENT: RF FULL - DEMOGRAPHY, PREDICITIONS
tstamp <- "20231002191703ywh7g"

predicted <- fread(paste0("output/predict/",tstamp,".predicted.csv"),data.table=F)


predicted_n <- data.frame(
    data=convertDemoStrings(predicted$data)[,"demo_3class_NcSi"],
    prediction=convertDemoStrings(predicted$prediction)[,"demo_3class_NcSi"]) %>%
    mutate_at(c("prediction","data"),as.factor)

cm <- caret::confusionMatrix(predicted_n$prediction,predicted_n$data)

f2s5 <- conf_matrix(cm=cm)
f2s5

## 2G: RF FULL - 3 CLASS DEMOGRAPHY, VARIABLE IMPORTANCE
tstamp <- "20231002191859ra2Tm"

varT <- var %>%
    filter(
        V9==tstamp,
        !is.na(imp)                                                                    
    )

f2s6 <- varplotFun(varT)
f2s6

## SUPPLEMENT: RF FULL - 3 CLASS DEMOGRAPHY, PREDICITIONS
tstamp <- "20231002191859ra2Tm"
tstampNew <- "2023100316274058c2c538"

predicted <- fread(
    paste0("output/predict/",tstamp,".run",tstampNew,".predicted.csv"),
    data.table=F)

predicted_n <- data.frame(
    data=demo3Function(predicted$data),
    prediction=demo3Function(predicted$prediction)) %>%
    mutate_at(c("prediction","data"),as.factor)

cm <- caret::confusionMatrix(predicted_n$prediction,predicted_n$data)

f2s7 <- conf_matrix(cm=cm)
f2s7


## SUPPLEMENT FACETED VAR IMP
varT <- var %>%
    filter(
        V9!="202310021917048AtJy", # sequential model
        !is.na(imp)                                                                    
    )

f2s8 <- varT %>%
    mutate(response=ifelse(V3=="selfRate","Selfing Rate",
                    ifelse(V3=="demo_3class","Demography\nbroad category",
                    ifelse(V3=="demography","Demography",NA)))) %>%
    varplotFun()+
    facet_grid(response~model,#switch="y",
               scales="free",space="free"
               )
f2s8

## SUPPLEMENT RF-FULL PREDICITIONS with selection visible
## 2E: RF FULL TEST-UNSEEN, PREDICITIONS
tstamp <- "20231002191703IfZEj"
tstampNew <- "202310031627384c9524cd"

predicted <- fread(
    paste0("output/predict/",tstamp,".run",tstampNew,".predicted.csv"),
    data.table=F)

load(paste0("output/predict/",tstamp,".run",tstampNew,".dTest.RData"))

predictedParams <- bind_cols(dTest,predicted) %>% # add parameters to prediction
    select(data,prediction,neutral) %>%
    mutate(Mutations=ifelse(neutral,"neutral","with selection"))

f2s9 <- ggplot(predictedParams)+ # from the viol_plotter function, with added parameters
    geom_violin(aes(x=paste(data),y=prediction,fill=Mutations))+
    geom_hline(aes(yintercept=data),color='red')+
    guides(x=guide_axis(angle=45))+
    labs(x="True Selfing",y="Prediction")+
    facet_grid(.~paste(data),scales="free",space="free")+
    theme(strip.text.x = element_blank(),
          panel.spacing = unit(0,"lines"))+
    background_grid(major="y")
f2s9

## Assemble final Figure 2
tp <- c(-.05,1)
f2b <-  f2b+ theme(plot.tag.position = tp)
f2d <-  f2d+ theme(plot.tag.position = tp)
f2s6 <- f2s6+theme(plot.tag.position = tp)

top <-  f2b|  # ROH-model variable importance
        f2c   # ROH-model unseen predictions
mid <-  f2d|  # full-model variable importance
        f2e   # full-model unseen predictions
bot <-  f2s6| # full-model 3 class demography variable importance
        f2s7  # full-model 3 class demography prediction

f2_plots <- top/mid/bot

title.properties <- ttheme_default(
    base_size = 15,
    core = list(fg_params=list(rot=-90,
                               fontface="bold")))
top_text <- ggdraw(tableGrob("RF-ROH: Selfing Rate",theme=title.properties))
mid_text <- ggdraw(tableGrob("RF-full: Selfing Rate",theme=title.properties))
bot_text <- ggdraw(tableGrob("RF-full: Demography",theme=title.properties))
all_text <- top_text/mid_text/bot_text

f2 <- patchwork::free(p_workflow)|f2_plots

f2 <- f2 +
    plot_annotation(tag_levels = 'A')+
    plot_layout(widths = c(6,6),tag_level = 'new')
f2 <- plot_grid(f2,all_text,rel_widths=c(1,0.02),labels=c("A",NA))
f2

saveplot(f2,"plots/20240220_f2",12,8)

## save supplements
saveplot(f2s3,"plots/20240220_f2s3",6,4) # new
saveplot(f2s4,"plots/20231009_f2s4",6,4)
saveplot(f2s5,"plots/20231009_f2s5",10,6)
saveplot(f2s8,"plots/20240316_f2s8",10,6) # new 2
saveplot(f2s9,"plots/20231009_f2s9",6,4)



## FIGURE 3; PCA, Loadings, Scree
library(recipes)
library(tidymodels)

tstamp <- "20231002191703IfZEj" # REG, FULL

load(paste0("output/predict/",tstamp,".dTrain.RData"))

pcaData <- dTrain

pcaRec <- recipe(
    pcaData,
    formula=selfRate~lengthMedian+lengthVar+roh_count_ind+countVar+proportionInROH+propVar+gapMedian+gapVar+Fis+TajimaD+demo_3class,
    ) %>%
    update_role(demo_3class,new_role = "id") %>% 
    step_pca(all_predictors(), id = "pca", keep_original_cols = TRUE)

pcaPrep <- prep(pcaRec)
pca_loading <- tidy(pcaPrep, id="pca")
pca_variances <- tidy(pcaPrep, id = "pca", type = "variance")

pca_bake <- bake(pcaPrep,
                 pcaData
                 )

pca_plotbase <- pca_bake %>%
    mutate(demo_3class_nice=demo3Function(demo_3class)) %>% 
    ggplot(aes(
        PC1, PC2,
        color=selfRate
    )) +
    geom_point(aes(
        shape=demo_3class_nice
    ),
    alpha = 0.6, size = 2)+
    scale_color_viridis_c(
        option = "plasma",
        direction = -1,
        end = .95,
        breaks = c(0,.5,1))

## add loadings
pca_loadings <- tidy(pcaPrep, id = "pca", type = "coef")
pca_loadings <- pivot_wider(pca_loadings,names_from="component",values_from="value")

pca_loadings <- inner_join(pca_loadings,
                           data.frame(terms=names(varnamesShort),varnames=varnamesShort),
                           by='terms')

pcf <- 12
pca_plotLoadings <- pca_plotbase +
    geom_segment(
        data = pca_loadings,
        aes(x = 0, y = 0, xend = PC1 * pcf, yend = PC2 * pcf),
        arrow = arrow(length = unit(0.02, "npc")),
        color = "red"
    ) +
    geom_label_repel(
        data = pca_loadings,
        aes(x = PC1 * pcf, y = PC2 * pcf,
            label=varnames),
        color="black",
        point.padding = 0.25,
        fill=alpha("white",.7),
        parse=T
    )+
    coord_cartesian(
        ylim=c(-12,NA),
        xlim=c(-4,6.5)
    )+
    labs(shape="Demography",
         color="Selfing Rate   ")+
    theme(
        legend.position='bottom',#c(.8,.8)
        legend.spacing.x = unit(.01,"npc"),
        legend.key.width = unit(.05,"npc"),
        legend.box="vertical"
        )+
    guides(shape = guide_legend(override.aes = list(size=5)),
           color = guide_colourbar())
pca_plotLoadings

## Variance explained/scree plot
pca_variances <- tidy(pcaPrep, id = "pca", type = "variance")

scree_plt <- pca_variances %>%
    filter(terms == "percent variance") %>%
    ggplot(aes(component, value)) +
    geom_col(just=0.5,fill="#372F60") +
    labs(x = "Principal Components", y = "Variance explained (%)")+
    scale_x_continuous(breaks=1:10)+
    background_grid(major="y",minor="y")
scree_plt

## PCA with empirical data - projected onto the same axes
load_params <- function(path) {
    f <- fread(path,
               data.table=F) %>%
        filter(chr!="chr") %>%
        mutate_at(c("lengthVar", "countVar", "propVar", "roh_count_ind",  "proportionInROH",
                    "lengthMedian", "gapMedian", "gapVar", "TajimaD",  "Fis", "selfRate"),
                  as.numeric)
    return(f)
}
load("../models/rohSelfing_Full.RData")

## - Arabis
params_arabis_raw <- load_params(
    paste0("output/empirical/selfrates/arabis/params.bcfvcf_",tstamp,".txt")
)
params_arabis_raw$pop <- substr(params_arabis_raw[,ncol(params_arabis_raw)],1,2)
params_arabis_raw <- filter(params_arabis_raw,!pop%in%c("GE", "NU", "RI", "TV"))
params_arabis <- predict(model[[3]],params_arabis_raw) # preprocessing step from sims
pops <- params_arabis$pop
pca_arabis <- bake(pcaPrep,data.frame(params_arabis,demo_3class=pops)) # This means that the PCA recipe retains the parameters, such as loadings and mean centering, from the original PCA analysis.
pca_arabis$demo_3class <- pops

## - Lyrata
params_lyrata_raw <- load_params(
    paste0("output/empirical/selfrates/alyr/params.ragtag_",tstamp,".txt")
)
params_lyrata_raw <- params_lyrata_raw %>% filter(chr!="scaffold_96_RagTag")
params_lyrata_raw$pop <-
    gsub(paste0("selfrates_|","|_|",tstamp,"|.params|.txt"),"",
         params_lyrata_raw[,ncol(params_lyrata_raw)])
params_lyrata <- predict(model[[3]],params_lyrata_raw)
namesDF <- data.frame(pop=c("MW0079WO", "MW0079WS",  "MW015", "NT1", "NT12", "NT8"), # fix lyrata names
                      selfing=c("outcrossing","selfing","outcrossing","selfing","outcrossing","outcrossing"),
                      popNew = c("MW7O","MW7S","MW1","NT1", "NT12", "NT8"))
params_lyrata <- inner_join(params_lyrata,namesDF) %>% select(-pop,pop=popNew)
pops <- params_lyrata$pop
pca_lyrata <- bake(pcaPrep,data.frame(params_lyrata,demo_3class=pops))
pca_lyrata$demo_3class <- pops

## - Thaliana
params_thaliana_raw <- load_params(
    paste0("output/empirical/selfrates/thal/params.TAIR10_",tstamp,".txt")
)
params_thaliana_raw$pop <-
    gsub(paste0("selfrates_|","|_|",tstamp,"|.params|.txt"),"",
         params_thaliana_raw[,ncol(params_thaliana_raw)])
params_thaliana <- predict(model[[3]],params_thaliana_raw)
pops <- params_thaliana$pop
pca_thaliana <- bake(pcaPrep,data.frame(params_thaliana,demo_3class=pops))
pca_thaliana$demo_3class <- pops

## - Worms
params_worms_raw <- load_params(
    paste0("output/empirical/selfrates/worms/params.worms_",tstamp,".txt")
)
params_worms_raw$pop <-
    gsub(paste0("selfrates_|worms","|_|",tstamp,"|.params|.txt"),"",
         params_worms_raw[,ncol(params_worms_raw)])
params_worms <- predict(model[[3]],params_worms_raw)
pops <- params_worms$pop
pca_worms <- bake(pcaPrep,data.frame(params_worms,demo_3class=pops))
pca_worms$demo_3class <- pops

## - merge everything
pca_everything <- bind_rows(
    data.frame(pca_arabis,set="A. alpina"),
    data.frame(pca_lyrata,set="A. lyrata"),
    data.frame(pca_thaliana,set="A. thaliana"),
    data.frame(pca_worms,set="C. briggsae")
)

annotation <- pca_everything %>%
    group_by(set,demo_3class) %>%
    summarise(PC1=mean(PC1),
              PC2=mean(PC2))
pca_plotbase_everything <- pca_everything %>%
    ggplot(aes(
        PC1, PC2,
        color=set
    )) +
    geom_point(
        alpha = 1, size = 2)+
    geom_label_repel(
        data=annotation,
        aes(label=demo_3class),
        color="black",
        min.segment.length = Inf,
        max.overlaps = Inf,
        box.padding = .1,
        fill=alpha("white",.5),
        label.padding = .1,
        force=2
    )+
    scale_color_brewer(type="qual",palette=2)+
    theme(legend.position=c(.65,.25),
          legend.text = element_text(face='italic'))+
    labs(color="Species")
pca_plotbase_everything

## Assemble Figure 3
f3 <- plot_grid(
    pca_plotLoadings,
    plot_grid(
        scree_plt,
        pca_plotbase_everything,
        labels=c("B","C"),
        nrow=2,
        align='v'
    ),
    labels=c("A",NA),
    ncol=2
)
f3

saveplot(f3,"plots/20240220_f3",12,8)

#### FIGURE 4A Prop vs count - empirical data
require(ggtext)

ix <- "proportionInROH"
iy <- "roh_count_ind"

plot_features_empirical_data <- bind_rows(
    data.frame(params_lyrata_raw,set="A. lyrata") %>%
    inner_join(namesDF) %>% select(-pop,pop=popNew),
    data.frame(params_arabis_raw,set="A. alpina"),
    data.frame(params_thaliana_raw,set="A. thaliana"),
    data.frame(params_worms_raw,set="C. briggsae")
)

plot_features_empirical <- plot_features_empirical_data %>%
    group_by(pop,set) %>%
    summarise(ix = mean(!!sym(ix)),
              iy = mean(!!sym(iy))) %>% 
    ggplot(
        aes(
            ix,
            iy
        )
    )+
    background_grid()+
    geom_point(
        aes(
            ## shape=N_extra
            color=set
        ),
        alpha=.75,
        size=3
    )+
    geom_text_repel(
        aes(label=pop)
    )+
    scale_color_brewer(type="qual",palette=2)+
    labs(x=parse(text=eval(names(varnames)[varnames==ix])),
         y=names(varnames)[varnames==iy])+
    scale_y_log10()+
    scale_x_continuous(limits=c(0.1,1.05),breaks=seq(.25,1,.25))+
    annotation_logticks(sides = "l")+
    theme(legend.position='none')
plot_features_empirical

## FIGURE 4B; Sequential model VAR IMP, PREDICTIONS (N), EMPIRICAL PREDICTIONS

load("geos.RData") # init data for empirical data
load("geos_regionColors.RData")

modelIds <- c("20231002191703IfZEj", # full RF
              "202310021917048AtJy") # sequential RF

plist <- list() # list for empirical plots based on different modelIds
selfData_Models <- data.frame()

for (modelId in modelIds) {
    selfData <- fread(paste0("output/empirical/selfrates/arabis/bcfvcf_",modelId,"_selfrates.txt"),data.table=F)

    selfData$pop <- substr(selfData$V2,1,2)

    selfDataChr <- left_join(selfData,geos,by="pop")
    selfData_ArabisPops <- selfDataChr %>%
        group_by(pop,region,region2.order,country) %>%
        summarise(selfing_rate=mean(V1),
                  lon=mean(lon),
                  lat=mean(lat)) %>%
        filter(!region%in%c('Scandinavia')) %>%
        mutate(selfing=ifelse(country%in%c("Italy","Greece"),"outcrossing","selfing"))
    selfData <- selfData_ArabisPops %>%
        mutate(pop=region) %>%
        select(-country,-lon,-lat)

    selfData$popC <-
        factor(selfData$pop,unique(selfData$pop[order(selfData$region2.order)]),ordered=T)

    selfDataChr$popC <-
        factor(selfDataChr$pop,unique(selfDataChr$pop[order(selfDataChr$region2.order)]),ordered=T)

    selfData_Arabis <- selfData


    ## numbers
    selfData_Arabis %>% group_by(region) %>% summarise(mean(selfing_rate))

#### 3D - Boxplot Selfing Rates A. lyrata
    selfData <- fread(paste0("output/empirical/selfrates/alyr/ragtag_",modelId,"_selfrates.txt"),data.table=F)

    selfData <- selfData[1:(nrow(selfData)-1),] # extra scaffold we don't want (scaffold_96_RagTag)

    selfData$pop <- gsub(paste0("selfrates_|","|_|",modelId,"|.txt"),"",selfData$V2)
    selfData$selfing_rate <- as.numeric(selfData$V1)

    selfData <- inner_join(selfData,namesDF) %>% select(-pop,pop=popNew)

    selfData_Lyrata <- selfData

    ## numbers
    selfData_Lyrata %>% group_by(pop) %>% summarise(mean(selfing_rate))


    ## Thaliana
    selfData <- fread(paste0("output/empirical/selfrates/thal/TAIR10_",modelId,"_selfrates.txt"),data.table=F)

    selfData$pop <- gsub(paste0("selfrates_|","|_|",modelId,"|.txt"),"",selfData$V2)
    selfData$selfing_rate <- as.numeric(selfData$V1)
    selfData$selfing <- "selfing"

    selfData_thaliana <- selfData

    ## numbers
    selfData_thaliana %>% group_by(pop) %>% summarise(mean(selfing_rate))

    ## Worms
    selfData <- fread(paste0("output/empirical/selfrates/worms/worms_",modelId,"_selfrates.txt"),data.table=F)

    selfData$pop <- "<i>C. briggsae</i>"
    selfData$selfing_rate <- as.numeric(selfData$V1)
    selfData$selfing <- "selfing"

    selfData_worms <- selfData

    ## numbers
    selfData_worms %>% group_by(pop) %>% summarise(mean(selfing_rate))

    ## merge selfing data
    selfData_All <- bind_rows(
        selfData_Arabis %>% mutate(Species="A. alpina"),
        selfData_Lyrata %>% mutate(Species="A. lyrata"),
        selfData_thaliana %>% mutate(Species="A. thaliana"),
        selfData_worms%>% mutate(Species=""),
        ) %>%
        mutate(Species=factor(Species,
                              levels=c("A. alpina", "A. lyrata", "A. thaliana","")))

    selfData_Models <- rbind(selfData_Models,
                             selfData_All %>% mutate(modelId))

    plist[[modelId]] <- selfData_All %>%
        ggplot()+
        geom_boxplot(aes(x=pop,y=selfing_rate,fill=Species),alpha=.8)+
        stat_summary(aes(x=pop,y=selfing_rate,group=selfing,fill=Species),
                     fun="mean",geom="point",shape=23,size=3)+
        scale_fill_brewer(type="qual",palette=2)+
        facet_grid(~Species,space="free",scales="free",switch="x")+
        labs(y="Selfing Rate",x="Population")+
        background_grid(major="y")+
        theme(legend.position='none',
              strip.text = element_text(face = "italic"),
              axis.title.x=element_blank(),
              panel.spacing = unit(0.1, "lines"),
              axis.text.x=ggtext::element_markdown())+
        guides(x=guide_axis(angle=30))+
        coord_cartesian(ylim=c(0,1))
}

saveplot(plist[["20231002191703IfZEj"]],"plots/selfRatesEmpirical_fullRF")

bot <- plist[["202310021917048AtJy"]]

## Assemble figure 4
top <- plot_features_empirical

f4 <- top/bot
f4 <- f4 + plot_annotation(tag_levels = 'A')
f4

saveplot(f4,"plots/20240220_f4",8,8)


## SUPPLEMENTAL PLOT: Fis Model Empirical Data
predFis <- plot_features_empirical_data %>%
    left_join(geos,by="pop") %>%
    mutate(pop=ifelse(set=="A. alpina",as.character(region),pop)) %>% 
    select(Fis,set,selfing,pop) %>%
    mutate(selfing_rate=2*Fis/(1+Fis),
           Species=set) %>%
    unique()

p <- ggplot(predFis)+
        geom_boxplot(aes(x=pop,y=selfing_rate,fill=Species),alpha=.8)+
        stat_summary(aes(x=pop,y=selfing_rate,group=selfing,fill=Species),
                     fun="mean",geom="point",shape=23,size=3)+
        scale_fill_brewer(type="qual",palette=2)+
        facet_grid(~Species,space="free",scales="free",switch="x")+
        labs(y="Selfing Rate",x="Population")+
        background_grid(major="y")+
        theme(legend.position='none',
              strip.text = element_text(face = "italic"),
              axis.title.x=element_blank(),
              panel.spacing = unit(0.1, "lines"),
              axis.text.x=ggtext::element_markdown())+
        guides(x=guide_axis(angle=30))
p

saveplot(p,"plots/FisEmpiricalPreds")


