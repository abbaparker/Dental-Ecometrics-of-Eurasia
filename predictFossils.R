# Predict target environmental variable from your ensemble of saved models
# trained across cross-validation folds

predict_saved_model <- function(
    model,
    model_type,
    newdata
){
  
  if(model_type %in% c("GLM","GAM","SVM","MARS")){
    return(as.vector(
      predict(model, newdata = newdata)
    ))
  }
  
  if(model_type == "RF"){
    return(as.vector(
      predict(
        model,
        data = newdata
      )$predictions
    ))
  }
  
  if(model_type == "GBM"){
    return(as.vector(
      predict(
        model,
        newdata = newdata,
        n.trees = model$n.trees
      )
    ))
  }
  
}
