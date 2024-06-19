## rerun aggregateResults.R with these params and adjustments
simSet = "unseen"
statsOutname <- "unseen_summary.stats.20240423.wmixedpops.RData"
rohDataFileName <- "unseen_summary.20240423.roh.txt.RData"

##... do not filter out non-P1 twoPop demographies

stats <- statsAll %>%
    group_by(demography,neutral,mutationRate,mating,rep) %>%
    summarise(TajimaD=mean(TajimaD),
              Fis=mean(Fis),
              n=n()) %>% # min number of reps
    ungroup %>%
    filter(n>=20) %>% select(-n)

##... rest is unchanged


## CALC RMSE and PLOT DATA
## results could be presented with predictions known vs predicted >>by demography<<

require(recipes)
require(tidymodels)

## 202310021917048AtJy 202404241145030f97aa53
## 20231002191703bwG5E 2024042411450215851c60
## 20231002191703IfZEj 2024042411450287fc8226

tstamps <- c(
    selfRate.rf.ROHstat=   "20231002191703IfZEj",
    selfRate.rf.ROH=       "20231002191703bwG5E",
    selfRate.rf.sequential="202310021917048AtJy"
)

tstampsx <- c( # in the same order
    "2024042411450287fc8226",
    "2024042411450215851c60",
    "202404241145030f97aa53"
)

comparePopStructure <- lapply(seq_along(tstamps), function(x){
    load(paste0("workspaces/",tstamps[x],".run",tstampsx[x],".dTest.RData"))
    pred <- fread(paste0("workspaces/",tstamps[x],".run",tstampsx[x],".predicted.csv"),
                  data.table=F) %>%
        rename(pred.data=data,pred.prediction=prediction)
    t_base <- t <- bind_cols(dTest,pred)
    rmse <- rmse(data.frame(t=t$pred.data,
                            p=t$pred.prediction),
                 truth=t,
                 estimate=p)$.estimate
    t_all <- t %>% mutate(rmse=rmse,compare="entire dataset")
    t <- t_base %>% filter(grepl("twoPop",demography),grepl("P1",demography))
    rmse <- rmse(data.frame(t=t$pred.data,
                            p=t$pred.prediction),
                 truth=t,
                 estimate=p)$.estimate
    t_mixed <- t %>% mutate(rmse=rmse,compare="admix./isol.")
    t <- t_base %>% filter(grepl("twoPop",demography),!grepl("P1",demography))
    rmse <- rmse(data.frame(t=t$pred.data,
                            p=t$pred.prediction),
                 truth=t,
                 estimate=p)$.estimate
    t_onlystructure <- t %>% mutate(rmse=rmse,compare="substructured admix./isol.")
    df <- bind_rows(t_all %>% mutate(model=names(tstamps[x])),
                    t_mixed %>% mutate(model=names(tstamps[x])),
                    t_onlystructure %>% mutate(model=names(tstamps[x])))
    return(df)
})

save(comparePopStructure,
     file=paste0("caretmodels/evaldata/selfRate.comparePopStructure.20240423.RData"))

load("output/predict/selfRate.comparePopStructure.20240423.RData")


cps <- do.call(bind_rows,comparePopStructure)

nms <- c("selfRate.rf.ROH", "selfRate.rf.ROHstat", "selfRate.rf.sequential" )
names(nms) <- c("RF ROH","RF full","RF sequential")
cps$model <- forcats::fct_recode(cps$model,!!!nms)

nms <- c("entire dataset", "substructured admix./isol.", "admix./isol.")
cps$compare <- factor(cps$compare,levels=nms)
names(nms) <- c("entire dataset\n(substructured)", "sample metapopulation with substructure",
                "sample within subpopulation")
cps$compare <- forcats::fct_recode(cps$compare,!!!nms)
cps <- cps %>%     filter(model=="RF full",compare!="entire dataset\n(substructured)")

anno <- cps %>%
    group_by(compare,model) %>%
    summarise(r=unique(paste("RMSE =",round(rmse,3))))

pLinesFacetx <- cps %>%
    ggplot()+
    geom_smooth(aes(x=selfRate.orig,
                    y=pred.prediction),method="lm",se=F)+
    geom_point(aes(x=selfRate.orig,
                   y=pred.prediction),
               size=.1)+
    ## stat_summary(aes(x=selfRate.orig,
    ##                  y=pred.prediction),
    ##              geom="errorbar",
    ##              fun.data=mean_cl_normal, # Hmisc::mean_sdl is double standard dev
    ##              width=.05,
    ##              color="steelblue")+
    ## background_grid()+
    geom_abline(slope=1,intercept=0,color="red",linetype="dashed",size=1)+
    ggh4x::facet_grid2(model~compare,
                       axes="all",
                       remove_labels="all")+
    coord_cartesian(ylim=c(-.05,1.05))+
    theme(legend.position='none')+
    theme(axis.text.x = element_text(angle = 30, hjust = .8))+
    labs(x="True Selfing Rate",y="Predicted Selfing Rate")+
    geom_label(data=anno,aes(x=.3,y=1,label=r))
pLinesFacetx

saveplot(pLinesFacetx,paste0("plots/allModels.pLinesFacet.substructured.20240424"),10,6)

#paste(t$demography,t$subpop)%>%unique
