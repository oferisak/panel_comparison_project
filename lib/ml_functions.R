library(bonsai)

xgb_def <-
  function(train_tbl = data.frame(),
           grid_size = 10,
           hyperparam_opt = T) {
    model_definitions <- list()
    if (hyperparam_opt) {
      model_definitions[['spec']] <-
        boost_tree(
          trees = tune(),
          #tree_depth = tune(),
          min_n = tune(),
          loss_reduction = tune(),
          #sample_size=tune(),
          #      sample_size = tune(),# use this (reduced sample size to reduce overfitting)
          #mtry = tune(),
          #stop_iter = 20,
          learn_rate = tune()
        ) %>% set_engine("lightgbm") %>%  set_mode("classification")
      model_definitions[['parameters']] <-
        parameters(model_definitions[['spec']])
      #model_definitions[['parameters']] <-
      #  update(model_definitions[['parameters']], mtry = finalize(mtry(), train_tbl))
      model_definitions[['tune_grid_regular']] <-
        grid_regular(model_definitions[['parameters']], levels = 3)
      model_definitions[['tune_grid_latin_hypercube']] <-
        grid_latin_hypercube(model_definitions[['parameters']], size = grid_size)
    } else{
      model_definitions[['spec']] <-
        boost_tree() %>% set_engine("lightgbm") %>%  set_mode("classification")
      model_definitions[['parameters']] <-
        parameters(model_definitions[['spec']])
    }
    #model_definitions[['parameters']]<-update(model_definitions[['parameters']],learn_rate=learn_rate(range=c(-3,-1)),sample_size=sample_prop(range=c(0.1,0.7)),tree_depth=tree_depth(range=c(2,5)))
    return(model_definitions)
  }

glmnet_def<-function(train_tbl=data.frame(),grid_size=5,hyperparam_opt=T){
  model_definitions<-list()
  if (hyperparam_opt){
    model_definitions[['spec']]<- logistic_reg(penalty = tune(),mixture = 1)%>%set_engine("glmnet") %>% set_mode("classification")
    model_definitions[['parameters']]<-parameters(model_definitions[['spec']])
    model_definitions[['tune_grid_regular']] <- expand.grid(penalty=c(10^seq(-4,0, length.out = grid_size)))
    model_definitions[['tune_grid_latin_hypercube']] <-grid_latin_hypercube(model_definitions[['parameters']],size=grid_size)
  }else{
    model_definitions[['spec']]<-logistic_reg(penalty=0.0001,mixture = 1)%>%set_engine("glmnet") %>% set_mode("classification")
    model_definitions[['parameters']]<-parameters(model_definitions[['spec']])
  }
  return(model_definitions)
}

downsample_training_data<-function(train_data,outcome_column,under_ratio=1){
  ds_rec <- recipe(as.formula(glue('{outcome_column}~.')), data = train_data) %>%
    step_downsample(!!sym(outcome_column),under_ratio=under_ratio)
  p_re <- prep(ds_rec, retain=T) 
  train_downsample <- bake(p_re, new_data = NULL) 
  ncase_vs_ncontrols<-paste0(train_downsample%>%count(!!sym(outcome_column))%>%pull(n), collapse='/')
  message(glue('Downsampled training data.. now there are {nrow(train_downsample)} rows ({ncase_vs_ncontrols})'))
  return(train_downsample)
  
}

# This function trains the model using a group-wise cross validation, the user provides a group coulmn and 
# each unique value in that group is a cv sample
train_model<-function(train_data,model_def,formula,group='phenotype_name',grid_size=20){
  model_res<-list()
  model_res[['formula']]<-formula
  
  workflow_model <- workflow() %>%
    add_formula(form) %>%
    add_model(model_def$spec)
  
  #train_vfold <- vfold_cv(train_downsample, strata = has_pub,v = nfolds)
  train_vfold<-group_vfold_cv(train_data,group=group)
  message(glue('using {group} to split cv resamples - will result in {nrow(train_vfold)}-fold cross validation '))
  
  metrics_set<-metric_set(yardstick::pr_auc,yardstick::accuracy)
  library(doParallel)
  library(finetune)
  cl <- makePSOCKcluster(6)
  registerDoParallel(cl)
  
  tune_results <- workflow_model %>%
    tune_race_anova(resamples = train_vfold, #CV object
                    grid=grid_size,
                    param_info=model_def[['parameters']],
                    metrics = metrics_set,
                    control=control_race(save_pred = T,verbose = T,event_level='second'))
  
  stopCluster(cl)
  model_res[['tune_results']]<-tune_results
  model_res[['best_params']]<-select_best(tune_results, "pr_auc")
  model_res[['train_metrics']]<-collect_metrics(model_res[['tune_results']])
  final_wf <- finalize_workflow(
    workflow_model,
    model_res[['best_params']]
  )
  model_res[['final_wf']]<-final_wf
  model_res[['model_fit']]<-final_wf %>%
    fit(data = train_downsample)
  return(model_res)
}