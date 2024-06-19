## rerun aggregateResults10.R with these params and adjustments
## CALC RMSE and PLOT DATA
## results could be presented with predictions known vs predicted >>by demography<<

require(recipes)
require(tidymodels)

## 202310021917048AtJy 202404241145030f97aa53
## 20231002191703bwG5E 2024042411450215851c60
## 20231002191703IfZEj 2024042411450287fc8226

tstamps <- c(
    selfRate.rf.ROHstat=   "20231002191703IfZEj"
)

tstampsx <- c( # in the same order
    "20240426102909a32fdd56"
)

comparePopStructure <- lapply(seq_along(tstamps), function(x){
    load(paste0("output/predict/",tstamps[x],".run",tstampsx[x],".dTest.RData"))
    pred <- fread(paste0("output/predict/",tstamps[x],".run",tstampsx[x],".predicted.csv"),
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
    df <- bind_rows(t_all %>% mutate(model=names(tstamps[x])))
    return(df)
})


cps <- do.call(bind_rows,comparePopStructure)

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
    ## ggh4x::facet_grid2(model~compare,
    ##                    axes="all",
    ##                    remove_labels="all")+
    coord_cartesian(ylim=c(-.05,1.05))+
    theme(legend.position='none')+
    theme(axis.text.x = element_text(angle = 30, hjust = .8))+
    labs(x="True Selfing Rate",y="Predicted Selfing Rate")+
    geom_label(data=anno,aes(x=.3,y=1,label=r))
pLinesFacetx

saveplot(pLinesFacetx,paste0("plots/allModels.pLinesFacet.N10.20240425"),6,6)

#paste(t$demography,t$subpop)%>%unique


