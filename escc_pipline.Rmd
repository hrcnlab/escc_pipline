---
title: "Diagnostic model of escc pipline"
author: "hrcnlab"
date: "2022-11-21"
output: html_document
---

This code include the key parts of escc pipline in our paper, which contains three parts including model screening,feature selection and external diagnostic validation.

## Part1

### Machine learning for model determination

```{r}
#loading the packages
pacman::p_load(mlr,mlrMBO,tidyverse,ggplot2,mlrMBO)
```

```{r}
mlr_data_final %>% 
  janitor::clean_names()-> mlr_his2_clean
dim(mlr_his2_clean);str(mlr_his2_clean);table(mlr_his2_clean$his)
TransTask_raw <- makeClassifTask(data = mlr_his2_clean, 
                                 target = "his2")

# split data into train/test
split_desc = makeResampleDesc(method='Holdout', 
                              stratify = TRUE, split=0.75)
set.seed(2020)
split = makeResampleInstance(split_desc, task = TransTask_raw)
trainid = split$train.inds[[1]]; testid = split$test.inds[[1]]
train = mlr_his2_clean[trainid, ]; test = mlr_his2_clean[-trainid, ]
```

```{r}
# Define different Learners
learner_rpart <- makeLearner(
  id = "rpart",
  cl = "classif.rpart",
  predict.type = "prob"
)

learner_kknn <- makeLearner(
  id = "kknn",
  cl = "classif.kknn",
  predict.type = "prob"
)

# classic nbs from package e1071
learner_nbs <- makeLearner(
  id = "nbs",
  cl = "classif.naiveBayes",
  predict.type = "prob"
)

# Neural Network from nnet package
learner_nnet <- makeLearner(
  id = "nnet",
  cl = "classif.nnet",
  predict.type = "prob"
)
# randomFores from packages randomForest
learner_rf <- makeLearner(
  id = "rf",
  cl = "classif.randomForest",
  predict.type = "prob"
)
# logression
learner_glm <- makeLearner(
  id = "logre",
  cl = "classif.glmnet",
  predict.type = "prob"
)

learner_svm <- makeLearner(
  id = "svm",
  cl = "classif.svm",
  predict.type = "prob"
)

```

```{r}
define_lrn <- function(lrn) {
  lrn_define_filter <- makeFilterWrapper(
    learner = lrn,
    fw.method = "FSelectorRcpp_information.gain",
    more.args = list("FSelectorRcpp_information.gain" = list(
      equal = TRUE
    ))
  )

  inner <- makeResampleDesc("CV", iters = 5, stratify = T, predict = "test")

  define_ps <- makeParamSet(
    makeIntegerParam("fw.abs", lower = 5, upper = 25)
  )

  ctrl <- makeTuneControlGrid(resolution = 30)

  define_tuner <- makeTuneWrapper(lrn_define_filter,
    resampling = inner,
    par.set = define_ps,
    control = ctrl,
    measures = list(
      multiclass.au1u,
      logloss,
      ber, acc,
      bac, mmce,
      # multiclass.au1u,
      multiclass.aunu,
      # mmcegroupMean,
      timetrain
    ),
    show.info = TRUE
  )
}
```

```{r}
lrn_list <- list(
  learner_rpart, learner_glm, learner_kknn, learner_nbs,
  learner_nnet, learner_rf, learner_svm
)
lrn_list_define <- lapply(lrn_list, define_lrn)

res_instance <- makeResampleInstance("CV",
  iters = 10, # 10
  stratify = TRUE,
  predict = "both",
  task = TransTask_raw
)

library(parallelMap)
library(parallel)

set.seed(12345678)
parallelStartSocket(cpus = detectCores() - 1)
bmr <- mlr::benchmark(
  task = TransTask_raw,
  learners = lrn_list_define,
  resampling = res_instance,
  measures = list(
    logloss,
    acc,
    ber, bac, mmce,
    multiclass.au1u,
    multiclass.aunu,
    timetrain
  ),
  models = TRUE,
  show.info = TRUE
)
parallelStop()
```

```{r}
perf <- getBMRPerformances(bmr, as.df = TRUE)
p_logloss <- ggboxplot(perf,
  x = "learner.id",
  y = "logloss", color = "learner.id",
  add = c("means_se", "jitter"), order = lels_avg,
  add.params = list(color = "darkgray"),
  desc_stat = "mean_sd",
  error.plot = "errorbar",
  legend = "none"
) + rotate_x_text(angle = 45) +
  geom_hline(
    yintercept =
      mean(filter(perf, learner.id == "svm.filtered.tuned") %>% pull(logloss)), linetype = "longdash"
  ) +
  scale_x_discrete(labels = str_split(lels_avg, "\\.", simplify = T)[, 1])

p_logloss
```

## Part2

### Feature selection

```{r}
d# define funciton for lasso of feature selecting
LassoSub.bio<-
  function(k=1, Xdata, Ydata){
    set.seed(k)
    s=sample(nrow(Xdata), size=0.8*nrow(Xdata))
    Xsub=Xdata[s,]
    Ysub=Ydata[s]
    model.sub=cv.glmnet(x=Xsub, y=Ysub, 
                        alpha=1,
                        parallel = FALSE,
                        standardize = TRUE, # drop sclae
                        type.measure = "class",
                        family="binomial")
    coef.sub=coef(model.sub, 
                  s='lambda.min')
    beta <- names( coef.sub[apply(coef.sub != 0, 1, any),] )
    return(beta)
  }
```

```{r}
library(glmnet)
library(caret)
lasso.stab <- mclapply(1:100, 
                      FUN=LassoSub.bio,
                      Xdata=as.matrix(mlr_data_final[,-44]), 
                      Ydata=as.matrix(mlr_data_final[,44]) ,
                      mc.cores =ifelse(Sys.info()[['sysname']] == "Windows",
                                        yes = 1,
                                        no = (detectCores()-1) ) )
genes_frec<-enframe( table(unlist(lasso.stab)) ) # for bionomal
genes_frec<-genes_frec[order(-genes_frec$value),]
genes_frec<-genes_frec[-1,]
head(genes_frec)
```

Plot imporatance for fig5B

```{r}
#plot imporatance for fig

library(ggtext)
library(glue)
library(ggplot2)

highlight = function(x, pat, color="black", family="") {
  ifelse( x %in% pat , 
          glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

ggplot(genes_frec, 
       aes(x=reorder(name,value),
           y=value, 
           fill =reorder(value,-value) ) ) + 
  geom_col(show.legend = FALSE, alpha = 0.6 ) + 
  theme_bw() +
  scale_x_discrete( 
    labels= function(x)
    highlight(x,genes_jiao, "red") ) +
  theme_classic() + 
  theme(
    axis.text.y=element_markdown(margin = unit(c(0, 2, 0, 0), "mm"),
                                 angle = 0), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=10,face = "bold.italic" )+
  labs(title=" plot importance of genes",
       x ="Symbol", y = "Importance") +
  coord_flip()
```

## Part3

### Diagnostic model finding

```{r}
svm_Task <- makeClassifTask(
  data = mutil_auc_Data,
  target = "his2"
)
```

```{r}
split_desc <- makeResampleDesc(
  method = "Holdout",
  stratify = TRUE,
  split = 0.7
)
set.seed(202)
split <- makeResampleInstance(split_desc, task = svm_Task)
train <- split$train.inds[[1]]
test <- split$test.inds[[1]]
table(mutil_auc_Data[train, ]$his2)
```

```{r}
svm_ruc <- makeLearner(
    id = 'svm_mutli',
    cl = 'classif.svm',
    scale=FALSE, 
  predict.type = "prob")

kernels <- c("polynomial", "radial", "sigmoid")
svmParamSpace <- makeParamSet(
  makeDiscreteParam("kernel", values = kernels),
  makeIntegerParam("degree", lower = 1, upper = 3),
  makeNumericParam("cost", lower = 0.1, upper = 100),
  makeNumericParam("gamma", lower = 0.1, 10))


randSearch <- makeTuneControlRandom(maxit = 200)

cvForTuning <- makeResampleDesc("CV", iters=5,
                                predict = "both",
                                stratify = T)

parallelStartSocket(cpus = detectCores()-1)
set.seed(2018)
tunedSvmPars <- tuneParams(svm_ruc, 
                           task = subsetTask(svm_Task, train),
                           resampling = cvForTuning,
                           par.set = svmParamSpace,
                           measures = list( 
                             multiclass.au1u,
                              acc,ber,
                                bac,mmce,logloss,
                                multiclass.au1u,
                                multiclass.aunu,
                                timetrain),
                     control = randSearch)

parallelStop()
```

```{r}
tunedSvm <- setHyperPars(svm_ruc,
                         par.vals = tunedSvmPars$x)
tunedSvmModel <- mlr::train(tunedSvm, svm_Task,subset=train)
train_svm <- predict(tunedSvmModel, task = svm_Task,subset=train)
calculateConfusionMatrix(train_svm, relative = FALSE,
                         sums = FALSE, set = "both")
```

```{r}
res_svm <- predict(tunedSvmModel, task = svm_Task,subset=test)

calculateConfusionMatrix(res_svm, relative = FALSE,
                         sums = FALSE, set = "both")

res_svm_trs<-res_svm$data %>% 
  arrange(truth)  %>% mutate(id=row_number()) %>% 
  dplyr::select(id,truth,"prob.Normal","prob.low_grade","prob.maglin",response)

svm.table<-with(res_svm_trs,table(res_svm_trs$truth, res_svm_trs$response))
svm.table

caret::confusionMatrix(res_svm_trs$response,res_svm_trs$truth)
```

####Extra IHC validate

```{r}
library(readxl)

dim(ihc_new)
table(ihc_new$groups, useNA = "always")

ihc_news <- ihc_new %>%
  as_tibble() %>%
  dplyr::select(c(genes_five, groups)) %>%
  mutate(type = case_when(
    groups == "Normal" ~ "Normal",
    groups == "low" ~ "low_grade",
    groups == "high" ~ "maglin",
    groups == "Cancer" ~ "maglin"
  ))
ihc_news$type <- factor(ihc_news$type,
  levels = c("Normal", "low_grade", "maglin")
)

ihc_res_three <- ihc_news %>%
  select(c(genes_three, "type")) %>%
  drop_na()

```

```{r}
dim(ihc_res_three);table(ihc_res_three$type)
```

```{r}
dataX <- as.data.frame(ihc_res_three[, -4])

ihc.pred <- predict(tunedSvmModel, newdata = dataX)

ihc_results <- ihc.pred$data
ihc_results$truth <- ihc_res_three$type
ihc_svm_table <- with(ihc_results, table(ihc_results$truth, ihc_results$response))
ihc_svm_table
ihc_svm_table_m <- caret::confusionMatrix(ihc_results$response, ihc_results$truth)
ihc_svm_table_m
```
