# Demo, developed for Palaeoverse lecture Nov. 19, 2026. Author: Abigail K. Parker
# code adapted from collaborative work, acknowledgements: 
# Indrė Žliobaitė,
# Miikka Tallavaara
# Liping Liu
# Juha Saarinen
# Jeremias Glöggler
# Oscar Wilson
# Alexander Bakhia
# Sania Zubaid 

# Cite as: Parker et al. (with above authors, in prep.) "Robust trait-based inference of 
# environmental conditions through time" [to be updated on publication]

# All you need for ecometrics ####
# Libraries used
library(tidyverse)
library(rnaturalearth)
library(sf)
library(mgcv)
library(earth)
library(e1071)
library(kernlab)
library(randomForest)
library(gbm)
library(scico)
library(plotmo)
theme_set(theme_bw())

# Model training, modern data ####
# Load data  
setwd("")  # set to folder where you downloaded: DentalTraits_Genus.csv
# EcoISEA#H_Eurasia_climate.csv
# Europe_GenusOccurrences_PPPA.csv
# makeEcomDataGenus.R
# runCrossValidation.R
# predictFossils.R

# 1) Species occurrences
# let's use large mammals in taxa Perissodactyla, Primates, Proboscidea, Artiodactyla
# This is a matrix of their occurrences at sample points, filtered to a box around ~Europe
eur_occs <- read.csv("Europe_GenusOccurrences_PPPA.csv", row.names = 1)

# 2) Species traits -- Provided dental traits of large herbivores
# read more: https://helda.helsinki.fi/bitstreams/f34f2830-c05f-4c54-96cd-dc12fcc65e57/download
traits <- read.csv("DentalTraits_Genus.csv")

# 3) Environmental data -- sampled to match occurrence samples
# Source: Eco-ISEA3H database, includes varied climate, geographic, vegetation data
clim <- read.csv("EcoISEA3H_Eurasia_climate.csv", row.names = 1)

# For mapping
defaultCRS <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
eurasia <- ne_countries(scale="medium", continent = c('europe', 'asia', "africa", "oceania"), returnclass = "sf")

clim_spatial <- sf::st_as_sf(clim, coords = c("X", "Y"), 
                              crs=defaultCRS) 
ggplot()+ geom_sf(data=eurasia)+  
  geom_sf(data=clim_spatial, aes(color=Tree_Mean), shape=18, size=4)+ #demo will look at tree cover %
  coord_sf(xlim=c(-12, 50), ylim= c(20,70), expand=F)+
  scale_color_scico(palette="turku", direction = -1)

# Build an ecometric data sheet with trait averages per community sample
## here is the necessary function
source("./makeEcomDataGenus.R")
# requires path to occurrence matrix, path MamOccs, and path to genus traits, pathEcoTraits
ecometric_sheet <- makeEcomData(pathMamOccs = "./Europe_GenusOccurrences_PPPA.csv", pathEcoTraits= "./DentalTraits_Genus.csv")

# Merge these results with observed climate variables
EurEco <- merge(ecometric_sheet, clim, by="HID")
EurEco <- read.csv("Europe_EcometricSheet.csv", row.names = 1)

# What does a trait mean look like on the map?
ecom_spatial <- sf::st_as_sf(EurEco, coords = c("X", "Y"), 
                             crs=defaultCRS) 

ggplot()+ geom_sf(data=eurasia)+  
  geom_sf(data=ecom_spatial, aes(color=HYP), shape=18, size=4)+
  coord_sf(xlim=c(-12, 50), ylim= c(20,70), expand=F)+
  scale_color_continuous(low = "#B0E2FF", high = "#6B2C91")

# Now, to model relationships between community traits and environment ####
# let's use these dental traits: "HYP", "ALX", "OL", "SF", "BUN", "OT", and predict tree cover 

# Model forms, one by one
# Set up predictor variables and parameters
all_predictors <- c("HYP", "ALX", "OL", "LOP", "BUN")

predictor_sets <- list(  # we can customize which predictors are used for which model form
    GLM  = c("HYP","ALX"),
    GAM  = c("HYP","ALX"),
    SVM  = all_predictors,
    MARS = all_predictors,
    RF   = all_predictors,
    GBM  = all_predictors)

#Formulas to run
make_standard_formula <- function(response, predictors){
  as.formula( paste(response,"~",paste(predictors, collapse = " + ") ))
}

make_glm_formula <- function(response, predictors, degree){
  glm_terms <- paste0( "poly(", predictors, ", ", degree, ")" )
  as.formula(paste(response,"~", paste(glm_terms, collapse = " + ") ))
}

make_gam_formula <- function(response, predictors){
  gam_terms <- paste0("s(",predictors, ", k=5)")
  as.formula(paste(response,"~", paste(gam_terms, collapse = " + ") ))
}

# some default parameters
glm_degree <- 1 # degree for polynomial regression

mars_degree <- 2 #MARS model allows interaction between 2 predictors

rf_trees <- 1000 # for Random Forest model
rf_mtry <- 1 # one variable split per node in RF

gbm_trees <- 2000 # for Gradient Boosting Machines
gbm_depth <-  3
gbm_shrinkage <-  0.02

# set environmental target variable
target <- "Tree_Mean"

# Linear Regression: start here to visualize single-predictor relationships
glm_formula <- make_glm_formula(target, predictor_sets$GLM, glm_degree)

glm_fit <- glm(glm_formula, data = EurEco, family = gaussian()) #fit model on all data
ggplot(EurEco, aes(x=HYP, y=Tree_Mean)) + geom_point()+
                               geom_smooth(method="lm", se=F, color="dodgerblue")
ggplot(EurEco, aes(x=ALX, y=Tree_Mean)) + geom_point()+
  geom_smooth(method="lm", se=F, color="dodgerblue")

# General Additive Model: smooth-curve relationships
gam_formula <- make_gam_formula(target,predictor_sets$GAM )

gam_fit <- mgcv::gam(gam_formula, data = EurEco) #fit model on all data
# inspect smooth relationships
plotmo(gam_fit)

# Gradient Boosting Machines: builds decision trees sequentially
# to reduce residuals from the previous tree
gbm_formula <- make_standard_formula(target, predictor_sets$GBM)

gbm_fit <- gbm::gbm( #fit model on all data
  formula = gbm_formula,
  data = EurEco,
  distribution = "gaussian",
  n.trees = gbm_trees,
  interaction.depth = gbm_depth,
  shrinkage = gbm_shrinkage,
  verbose = FALSE
)

# Variable influence in nodes fit for this model
# tells you the importance of each variable for describing your target
summary(gbm_fit)

# MARS: Multivariate Adaptive Regression Splines
# fits piecewise linear splines across predictor space, rather than single line
mars_formula <- make_standard_formula(target, predictor_sets$MARS)

mars_fit <- earth::earth(mars_formula, data = EurEco,
  degree = mars_degree )

#inspect piecewise relationships
plotmo(mars_fit)

# SVM: Support Vector Machines
# fits hyperplane in predictor space for prediction of target
svm_formula <- make_standard_formula(target, predictor_sets$SVM)
svm_fit <- ksvm(svm_formula,data = EurEco) #fit model on all data

# Random Forest: start here to learn about variable importance ('feature selection')
# Builds ensemble of independent decision trees
rf_formula <- make_standard_formula(target, predictor_sets$RF) 

rf_fit <- randomForest(
  formula = rf_formula,
  data = EurEco,
  num.trees = rf_trees,
  mtry = rf_mtry,
  seed = 56)
# inspect responses
plotmo(rf_fit)

# Which of these model forms should we use? ####
# Assess their performance using spatial cross-validation

# Let's divide Europe into 1,000km width hexagonal blocks
scv <- cv_spatial(x = ecom_spatial, size=1000000, k=5)  #2900km hexagon size 
cv_plot(cv = scv, x = ecom_spatial, ncol=2) #plots how hexagons are assigned to CV folds

# Assign each sample community to one of these 5 blocks
# We will leave one aside for testing in each cross-validation round
EurEco$FoldID <- scv$folds_ids
# For further testing robustness of relationships across space,
# also add more iterations of spatial blocking

# Here is the cross-validation function, using the same models and parameters as above
source("./runCrossValidation.R")

# Create a directory to save the models trained per folds
dir.create("models_trees",recursive = TRUE,showWarnings = FALSE)
dir.create("test_results",recursive = TRUE,showWarnings = FALSE)

# Run this function to save cross-validated models 
CV_results <- run_5fold_cv(EurEco, fold_col = "FoldID") #prints fold progress

# Build ensemble -- average prediction from model forms for all test points
ensemble_results <- CV_results |>
  group_by(HID, Test_Fold) |>
  summarise(Observed = first(Observed), Prediction = mean(Prediction),
    .groups = "drop") |>
  mutate(Model = "Ensemble", Anomaly = Prediction - Observed)

all_results <- bind_rows(CV_results,ensemble_results)

# Compare individual model forms and ensemble cross-validation performance
# this can also be done for individual folds
performance <- all_results |>
  group_by(Model) |>
  summarise( #various statistics summarizing performance for models
    MAE  = mean(abs(Anomaly), na.rm = TRUE),
    RMSE = sqrt(mean(Anomaly^2, na.rm = TRUE)),
    Bias = mean(Anomaly, na.rm = TRUE),
    SD   = sd(Anomaly, na.rm = TRUE),
    R2   = cor(Observed, Prediction, use = "complete.obs")^2,
    N = n(), .groups = "drop" )

# Violin plot of model error, labelled with RMSE and bias (towards overprediction or 
    # underprediction)
comp_viol<- ggplot(all_results, aes(x = Model,y = Anomaly,fill = Model)) +
  geom_violin(width=0.85) +
  geom_boxplot(width = 0.15,outlier.shape = NA) +
  labs(x = "Model Form", y = "Prediction Error", fill = "Model") +
  ylim(-60, 60)+ #-25, 25 for temp
  theme(legend.position = "none",
         axis.text.x = element_text(angle = 40,hjust = 1 ), color="black", size=11)+
  scale_fill_brewer(palette = "Set3",labels = c("GLM", "GAM", "MARS", "SVM", "RF", "GBM", "Ensemble")) +
  scale_x_discrete(limits = c("GLM", "GAM", "MARS", "SVM", "RF", "GBM", "Ensemble"))+
  geom_label(data = performance, aes(x = Model,y = -50,
                                     label = paste0(
                                       "Bias: ", format(signif(Bias, 2), scientific = TRUE),
                                       "\nRMSE: ", round(RMSE, 2) ) ),
             size = 3.2, linewidth = 0.1)

comp_viol #plot

# Save predictions across test folds
date_stamp <- format(Sys.Date(), "%Y-%m-%d")
write.csv(CV_results, file.path("./test_results/", paste0(target, "_CV_predictions_",
  date_stamp, ".csv")), row.names = FALSE)

# Apply to some fossil sites ####
# Here, trait averages have been calculated ahead of time for 
now_europe <- read.csv("NOWDatabase_Pleistocene_europe_sites.csv", row.names = 1)

# Locate models
model_files <- list.files("models_trees",
  pattern = "\\.rds$",recursive = TRUE,full.names = TRUE)
model_files <- model_files[!grepl("_SVM\\.rds$", model_files)] #here, we exclude SVM

# Function to predict from your saved models
source("./predictFossils.R")

prediction_list <- vector("list",length(model_files))

# Get predictions from 6x5 models
for(i in seq_along(model_files)){
  model_file <- model_files[i]
  cat("Predicting with", basename(model_file), "\n")
  
  model_name <- tools::file_path_sans_ext(basename(model_file))
  model_type <- sub(".*_","",model_name)
  
  test_fold <- as.numeric(sub(".*_Fold([0-9]+)_.*","\\1", model_name))
  fit <- readRDS(model_file)

  preds <- predict_saved_model(model = fit, model_type = model_type, newdata = now_europe)
  
  prediction_list[[i]] <- data.frame(
    HID = now_europe$HID,
    Test_Fold = test_fold,
    Model = model_type,
    Prediction = preds
  )
}
fossil_predictions <- bind_rows(prediction_list)

# Build ensemble prediction for each site
ensemble_predictions <- fossil_predictions %>%
  group_by(HID,Test_Fold) %>%
  summarise(Prediction = mean(Prediction),.groups = "drop") %>%
  mutate(Model = "Ensemble")

all_fossil_predictions <- bind_rows(
  fossil_predictions, ensemble_predictions)

# Save predictions
write.csv(all_fossil_predictions, (paste0( target,"_Predictions.csv")), 
  row.names  = FALSE)

# Summarize by site
fossil_summary <- all_fossil_predictions %>% 
  group_by(HID) %>%
  summarise(
    Mean = mean(Prediction, na.rm=T),
    Median = median(Prediction, na.rm=T),
    SD = sd(Prediction, na.rm=T),
    Lower95 = quantile(Prediction,0.025, na.rm=T),
    Upper95 = quantile( Prediction,0.975, na.rm=T),
    N = n(), .groups = "drop")

# Plot the range of reconstructions per site
# Find the overall mean across model forms and folds
# Small subset for plot
sites_sub <- fossil_summary %>% dplyr::slice_sample(n=25)
  
#Plot these tree percentage reconstructions
ggplot(sites_sub,
  aes(x = reorder(HID, Mean), y = Mean)) +
  geom_point(size=4, color="dodgerblue") +
  
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95),
    width = 0.3) +
  
  coord_flip() +
  
  theme_bw() +
  
  labs(
    x = "Fossil Site",
    y = target,
    title = paste(
      target,
      "reconstructions from models with spatial CV"
    )
  )



