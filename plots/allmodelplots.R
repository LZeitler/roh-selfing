#### THIS SCRIPT IS TO MAKE "allModels*" PLOTS, INCLUDING LASSO
#### FOR SUBSAMPLING AND POPSTRUCTRE RESULTS, RUN n10eval.R or admixpopstructure.R 
## to save time you can load an intermediate output file, search for "JumpHere" skip the time
## intensive stuff
## (this file is ex lasso.R)

source('source_me_light.R')

require(recipes)
require(tidymodels)
require(patchwork)

tstamp <- "20231002191703IfZEj" # REG, have to load full model here, otherwise not preprocessed
load(paste0("workspaces/workspace",tstamp,".RData"))


#### LASSO GLMNET MODEL - with penalty optimization using bootstrap

## prep recipe, dTrain is already preprocessed, i.e., normalized (center, scale)
lassoRec <- recipe(
    data=dTrain, 
    formula=as.formula(formulaStr)
)

## Specify model
lassoSpecGlmnet <- linear_reg(
    penalty=1,
    mixture=0.1
) %>%
    set_engine("glmnet")

## apply model to prepped data
wf <- workflow() %>%
    add_recipe(lassoRec)
lassoFitGlmnet <- wf %>%
    add_model(lassoSpecGlmnet) %>%
    fit(data=dTrain)

## look at results
lassoFitGlmnet %>%
    extract_fit_parsnip() %>% 
    tidy()

## RMSE
predictions <- lassoFitGlmnet %>% predict(dTrain)

rmse(tibble(t=dTrain[,responseVariable],p=predictions$.pred),
     truth=t,estimate=p)$.estimate



## tune penalty to optimize model
boots <- bootstraps(dTrain,strata=rep)

lassoSpecGlmnetTune <- linear_reg(
    penalty = tune(),
    mixture = 1
) %>%
    set_engine("glmnet")

lambda_grid <- grid_regular(penalty(), levels = 20)

lasso_grid <- tune_grid(
    wf %>% add_model(lassoSpecGlmnetTune),
    resamples = boots,
    grid = lambda_grid
)

## overview boot results
lasso_grid %>%
    collect_metrics()

## select best parameters and use for final model
good_rmse <- lasso_grid %>%
    select_best("rmse")

lassoSpecGlmnetFinal <- finalize_workflow(
  wf %>% add_model(lassoSpecGlmnetTune),
  good_rmse
)

## estimates of predictors, using final model
lassoSpecGlmnetFinal %>%
    fit(dTrain) %>%
    extract_fit_parsnip() %>% 
    tidy()

## apply model to dTrain
lfit <- lassoSpecGlmnetFinal %>%
    fit(dTrain)
predictions <- lfit %>% predict(dTrain)

rmse(tibble(t=dTrain[,responseVariable],p=predictions$.pred),
     truth=t,estimate=p)$.estimate

## apply model to dTest
predictions <- lfit %>% predict(dTest)

rmse(tibble(t=dTest[,responseVariable],p=predictions$.pred),
     truth=t,estimate=p)$.estimate


#### LINEAR REGRESSION
## Specify model
lregSpec <- linear_reg() %>% set_engine("lm")

## Prep recipes
lregRec1 <- recipe(
    data=dTrain, 
    formula=as.formula(paste0(responseVariable,"~Fis"))
)

lregRec2 <- recipe(
    data=dTrain, 
    formula=as.formula(paste0(responseVariable,"~lengthVar+proportionInROH"))
)

wfset <- workflow_set(
    preproc = list(Fis = lregRec1,
                   ROH = lregRec2), 
    models = list(lregSpec)
)

## Prep workflow, include recipes
grid_ctrl <- control_grid(
    save_pred = TRUE,
    save_workflow = TRUE
)

grid_results <- wfset %>%
    workflow_map(
        resamples = bootstraps(dTrain,0,apparent=T),
        control = grid_ctrl
    ) %>%
    mutate(
        fit = map(result, fit_best, metric = "rmse")
    )

## Looks at results
grid_results %>% 
    rank_results() %>% 
    filter(.metric == "rmse") %>% 
    select(wflow_id, model, .config, rmse = mean, rank)

## Fis_linear_reg
grid_results$result[[1]]$.predictions # prediction train set
lregFisPred <- grid_results$fit[[1]] %>% predict(dTest) %>%
    pull(.pred) # predictions test set

## ROH_linear_reg
grid_results$result[[2]]$.predictions
lregROHPred <- grid_results$fit[[2]] %>% predict(dTest) %>%
    pull(.pred)



#### POPULATION GENETIC FORMULA - uses Fis
#### S = 2F / (1 + F)
rmse(tibble(t=dTrain[,responseVariable],p=2*dTrain$Fis.orig/(1+dTrain$Fis.orig)),
     truth=t,estimate=p)
rmse(tibble(t=dTest[,responseVariable],p=2*dTest$Fis.orig/(1+dTest$Fis.orig)),
     truth=t,estimate=p)

#### APPLY MODELS TO DATA
#### Assemble metrics for dTest for all models
## linear models from above
allPredictions <- dTest %>%
    mutate(selfRate.Popgen=2*dTest$Fis.orig/(1+dTest$Fis.orig),
           selfRate.lreg.Fis=lregFisPred,
           selfRate.glmnet.ROHstat=predictions$.pred
           )

## rf models. The time stamps refer to models - the data has a constant tstamp, namely from the full model
tstamps <- c(
    selfRate.rf.ROHstat=   "20231002191703IfZEj",
    selfRate.rf.ROH=       "20231002191703bwG5E",
    selfRate.rf.stat=      "20231002191703RXUjU",
    selfRate.rf.sequential="202310021917048AtJy"
)
rf_results <- sapply(tstamps, function(x){
    load(paste0("caretmodels/model",x,".RData"))
    model[[1]] %>% predict(dTest)
})
allPredictions <- bind_cols(allPredictions,rf_results)

## save(allPredictions, file = paste0("caretmodels/evaldata/selfRate.allPredictions.",tstamp,".RData"))

#### JumpHere to save some time and memory
##--- LOCAL ---
load(paste0("output/selfRate.allPredictions.",tstamp,".RData"))

newPredictions <- allPredictions%>%select(contains("selfRate."),-selfRate.orig)

rmses <- sapply(1:ncol(newPredictions), function(x) {
    rmse(data.frame(t=allPredictions$selfRate.orig,
                    p=newPredictions[,x]),
         truth=t,
         estimate=p)$.estimate
})
names(rmses) <- names(newPredictions)
sort(rmses)

qplot(allPredictions$selfRate,allPredictions$selfRate.rf.ROHstat)

## plots
allPredLong <- allPredictions %>%
    select(contains("selfRate.")) %>%
    pivot_longer(!selfRate.orig,
                 names_to="model",
                 values_to="prediction",
                 names_prefix="selfRate.")

nms <- c("F[IS]~based"="Popgen", # for presentations
         "Linear~Regression~on~F[IS]"="lreg.Fis",
         "LASSO-GLMNET"="glmnet.ROHstat",
         "RF-sequential"="rf.sequential",
         "RF-full"="rf.ROHstat",
         "RF-ROH"="rf.ROH",
         "RF-stats"="rf.stat")


## reorder models
allPredLong$model2 <- factor(allPredLong$model,
                             levels=gsub("selfRate.","",names(sort(rmses))))

## rename models
allPredLong$model2 <- forcats::fct_recode(allPredLong$model2,!!!nms)

anno <- data.frame(r=rmses)
anno$model <- gsub("selfRate.","",row.names(anno))
anno$r <- paste("RMSE =",round(anno$r,3))
anno$model2 <- forcats::fct_recode(anno$model,!!!nms)

leg <- get_legend(data.frame(clr=c("steelblue","red"),
                             expl=c("linear fit","identity line")) %>%
                  ggplot()+
                  geom_histogram(aes(x=clr,fill=expl),stat="count")+
                  scale_fill_manual(values=rev(c("steelblue","red")))+
                  theme(legend.title = element_blank(),
                        legend.text = element_text(size=16)))

pLinesFacet1 <- allPredLong %>% 
    ggplot()+
    geom_smooth(aes(x=selfRate.orig,
                    y=prediction),method="lm",se=F)+
    geom_point(aes(x=selfRate.orig,
                   y=prediction),
               size=.1)+
    geom_abline(slope=1,intercept=0,color="red",linetype="dashed",size=1)+
    facet_wrap(.~model2,nrow=2, labeller = "label_parsed")+
    coord_cartesian(ylim=c(-.05,1.05))+
    theme(legend.position='none')+
    theme(axis.text.x = element_text(angle = 30, hjust = .8))+
    labs(x="True Selfing Rate",y="Predicted Selfing Rate")+
    geom_label(data=anno,aes(x=.3,y=1,label=r))
pLinesFacet1

pleg <- pLinesFacet1+leg+
    plot_layout(design=c(
                    area(1,1,20,20),
                    area(11,17,20,20)))
pleg

saveplot(pleg,paste0("plots/allModels.pLinesFacet.newnames.pres106",tstamp))


#### Apply models to unseen data
## Again you can jump to JumpHere
rohDataSummary <- stats <- NULL
statsStr <- "unseen_summary.stats.20231004byPopOnlyP1.RData"
load(statsStr)

rohStr <- "unseen_summary.20231004.roh.txt.RData"
load(rohStr)

dUnseen <- inner_join(rohDataSummary,stats) %>%
    mutate(selfRate=as.numeric(as.character(
               factor(mating,
                      levels=c("outc100","self100","self25","self75"),
                      labels=c(0,1,.25,.75))))) %>%
    ungroup

dUnseenTransformed <- predict(preProcValues, dUnseen)   # preproc from dTrain!
names(dUnseen) <- paste0(names(dUnseen), ".orig")
dUnseen <- bind_cols(dUnseen,
                     dUnseenTransformed)


## LASSO
pred_lasso <- lfit %>% predict(dUnseen) %>%
    pull(.pred)

## LReg Fis
pred_lregFis <- grid_results$fit[[1]] %>% predict(dUnseen) %>%
    pull(.pred)

## Popgen
pred_popgen <- 2*dUnseen$Fis.orig/(1+dUnseen$Fis.orig)

## everything
allPredictionsU <- dUnseen %>%
    bind_cols(
        selfRate.Popgen=         pred_popgen,
        ## selfRate.lreg.ROH=       pred_lregROH,
        selfRate.lreg.Fis=       pred_lregFis,
        selfRate.glmnet.ROHstat= pred_lasso,
    )

## RF models
newtstamps <- c(
    selfRate.rf.ROHstat="202310031627384c9524cd",
    selfRate.rf.ROH="202310031627367d95b179",
    selfRate.rf.stat="20231003162713d447a05c",
    selfRate.rf.sequential="202310041738061ea7edd7"
)
pred_rf <- lapply(names(tstamps), function(x){
    fread(
        paste0("workspaces/",tstamps[x],".run",newtstamps[x],".predicted.csv"),
        data.table=F)[,2]
})
names(pred_rf) <- names(tstamps)
pred_rf <- bind_cols(pred_rf)
allPredictionsU <- bind_cols(allPredictionsU,pred_rf) %>% data.frame

## save(allPredictionsU, file =
## paste0("caretmodels/evaldata/selfRate.allPredictionsU.",tstamp,".RData"))
## JumpHere for unseen data
load(paste0("output/selfRate.allPredictionsU.",tstamp,".RData"))


newPredictions <- allPredictionsU%>%select(contains("selfRate."),-selfRate.orig)

rmses <- sapply(1:ncol(newPredictions), function(x) {
    rmse(data.frame(t=allPredictionsU$selfRate.orig,
                    p=newPredictions[,x]),
         truth=t,
         estimate=p)$.estimate
})
names(rmses) <- names(newPredictions)
sort(rmses)

allPredLong <- allPredictionsU %>%
    select(contains("selfRate.")) %>%
    pivot_longer(!selfRate.orig,
                 names_to="model",
                 values_to="prediction",
                 names_prefix="selfRate.")

## reorder models
allPredLong$model2 <- factor(allPredLong$model,
                             levels=gsub("selfRate.","",names(sort(rmses))))

## rename models
allPredLong$model2 <- forcats::fct_recode(allPredLong$model2,!!!nms)

anno <- data.frame(r=rmses)
anno$model <- gsub("selfRate.","",row.names(anno))
anno$r <- paste("RMSE =",round(anno$r,3))
anno$model2 <- forcats::fct_recode(anno$model,!!!nms)

leg <- get_legend(data.frame(clr=c("steelblue","red"),
                             expl=c("linear fit","identity line")) %>%
                  ggplot()+
                  geom_histogram(aes(x=clr,fill=expl),stat="count")+
                  scale_fill_manual(values=rev(c("steelblue","red")))+
                  theme(legend.title = element_blank(),
                        legend.text = element_text(size=16)))

pLinesFacet1 <- allPredLong %>% 
    ggplot()+
    geom_smooth(aes(x=selfRate.orig,
                    y=prediction),method="lm",se=F)+
    geom_point(aes(x=selfRate.orig,
                   y=prediction),
               size=.1)+
    geom_abline(slope=1,intercept=0,color="red",linetype="dashed",size=1)+
    facet_wrap(.~model2,nrow=2, labeller = "label_parsed")+
    coord_cartesian(ylim=c(-.05,1.05))+
    theme(legend.position='none')+
    theme(axis.text.x = element_text(angle = 30, hjust = .8))+
    labs(x="True Selfing Rate",y="Predicted Selfing Rate")+
    geom_label(data=anno,aes(x=.3,y=1,label=r))
pLinesFacet1

pleg <- pLinesFacet1+leg+
    plot_layout(design=c(
                    area(1,1,20,20),
                    area(11,17,20,20)))
pleg

saveplot(pleg,paste0("plots/allModels.pLinesFacet.Unseen.newnames1.pres106",tstamp))
