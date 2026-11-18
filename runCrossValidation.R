# Using defined spatial folds, run 5-fold cross validation for our models

run_5fold_cv <- function(EurEco, fold_col = "FoldID"){
  
  all_predictions <- list()
  
  model_dir <- "models_trees"
  
  for(test_fold in 1:5){
    
    cat("Fold", test_fold, "\n")
    
    train_df <- EurEco[
      !is.na(EurEco[[fold_col]]) &
        EurEco[[fold_col]] != test_fold,
    ]
    
    test_df <- EurEco[
      !is.na(EurEco[[fold_col]]) &
        EurEco[[fold_col]] == test_fold,
    ]
    
    prediction_list <- list()
    
    # GLM --------------------------------------------------
    
    glm_formula <- make_glm_formula(
      target,
      predictor_sets$GLM,
      glm_degree
    )
    
    glm_fit <- glm(
      glm_formula,
      data = train_df,
      family = gaussian()
    )
    
    saveRDS(
      glm_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_GLM.rds")
      )
    )
    
    glm_pred <- predict(
      glm_fit,
      newdata = test_df
    )
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "GLM",
      Observed = test_df[[target]],
      Prediction = glm_pred,
      Anomaly = glm_pred - test_df[[target]]
    )
    
    rm(glm_fit)
    gc()
    
    # GAM --------------------------------------------------
    
    gam_formula <- make_gam_formula(
      target,
      predictor_sets$GAM
    )
    
    gam_fit <- mgcv::gam(
      gam_formula,
      data = train_df,
      method = "REML"
    )
    
    saveRDS(
      gam_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_GAM.rds")
      )
    )
    
    gam_pred <- predict(
      gam_fit,
      newdata = test_df
    )
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "GAM",
      Observed = test_df[[target]],
      Prediction = gam_pred,
      Anomaly = gam_pred - test_df[[target]]
    )
    
    rm(gam_fit)
    gc()
    
    # SVM --------------------------------------------------
    
    svm_formula <- make_standard_formula(
      target,
      predictor_sets$SVM
    )
    
    svm_fit <- ksvm(
      svm_formula,
      data = train_df
    )
    
    saveRDS(
      svm_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_SVM.rds")
      )
    )
    
    svm_pred <- predict(
      svm_fit,
      newdata = test_df
    )
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "SVM",
      Observed = test_df[[target]],
      Prediction = svm_pred,
      Anomaly = svm_pred - test_df[[target]]
    )
    
    rm(svm_fit)
    gc()
    
    # MARS -------------------------------------------------
    
    mars_formula <- make_standard_formula(
      target,
      predictor_sets$MARS
    )
    
    mars_fit <- earth::earth(
      mars_formula,
      data = train_df,
      degree = mars_degree
    )
    
    saveRDS(
      mars_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_MARS.rds")
      )
    )
    
    mars_pred <- predict(
      mars_fit,
      newdata = test_df
    )
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "MARS",
      Observed = test_df[[target]],
      Prediction = as.vector(mars_pred),
      Anomaly = as.vector(mars_pred) - test_df[[target]]
    )
    
    rm(mars_fit)
    gc()
    
    # RF ---------------------------------------------------
    
    rf_formula <- make_standard_formula(
      target,
      predictor_sets$RF
    )
    
    rf_fit <- ranger::ranger(
      formula = rf_formula,
      data = train_df,
      num.trees = rf_trees,
      mtry = rf_mtry,
      seed = 123
    )
    
    saveRDS(
      rf_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_RF.rds")
      )
    )
    
    rf_pred <- predict(
      rf_fit,
      data = test_df
    )$predictions
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "RF",
      Observed = test_df[[target]],
      Prediction = rf_pred,
      Anomaly = rf_pred - test_df[[target]]
    )
    
    rm(rf_fit)
    gc()
    
    # GBM --------------------------------------------------
    
    gbm_formula <- make_standard_formula(
      target,
      predictor_sets$GBM
    )
    
    gbm_fit <- gbm::gbm(
      formula = gbm_formula,
      data = train_df,
      distribution = "gaussian",
      n.trees = gbm_trees,
      interaction.depth = gbm_depth,
      shrinkage = gbm_shrinkage,
      verbose = FALSE
    )
    
    saveRDS(
      gbm_fit,
      file.path(
        model_dir,
        paste0(target, "_Fold", test_fold, "_GBM.rds")
      )
    )
    
    gbm_pred <- predict(
      gbm_fit,
      newdata = test_df,
      n.trees = gbm_trees
    )
    
    prediction_list[[length(prediction_list)+1]] <- data.frame(
      HID = test_df$HID,
      Test_Fold = test_fold,
      Model = "GBM",
      Observed = test_df[[target]],
      Prediction = gbm_pred,
      Anomaly = gbm_pred - test_df[[target]]
    )
    
    rm(gbm_fit)
    gc()
    
    all_predictions[[test_fold]] <-
      dplyr::bind_rows(prediction_list)
  }
  
  dplyr::bind_rows(all_predictions)
}