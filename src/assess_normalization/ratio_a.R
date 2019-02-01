# Robin Gradin
#
# Batch adjustment by reference alignment (2019)
# Prediction results for ratio-A (bapred)

# Load dependencies
library(bapred)
library(ggplot2)
library(ggthemes)
library(randomForest)
library(reshape2)
library(RGvisualization)
library(stringr)
library(viridis)

source('R code/bara.R')
source('R code/globals.R')
source('R code/helpers.R')
source('R code/performance_functions.R')
source('R code/prediction_fit_functions.R')
source('R code/selection_functions.R')
source('R code/optimize_parameters.R')
source('R code/assess_normalization.R')

# Load data
datasets <- readRDS(file = 'processed_datasets.RDS')

# Define parameters for the function
number_of_genes <- 500
fs_fun <- function(...){
  res <- select_significant_limma(...)
  return(res$selected_genes)
}

normalization <- list(
  fit = function(x, ref, ...){
    batch <- factor(rep(1, nrow(x)))
    fit <- bapred::ratioa(x = x, batch = batch)
    return(list(
      x_train = fit$xadj,
      fit = fit,
      ref = ref
    ))
  },
  transform = function(fit, x_test, ref){
    batch <- factor(rep(1, nrow(x_test)))
    res <- bapred::ratioaaddon(params = fit$fit, x = x_test, batch = batch)
    return(list(x_test = res, ref = ref))
  }
)

function_params_svm <- list(
  'fs_fun' = list(n = number_of_genes),
  'fit_fun' = list(kernel = 'linear'),
  'pred_fun' = list(decision.values = FALSE)
)
function_params_rf <- list(
  'fs_fun' = list(n = number_of_genes)
)

function_params_knn <- list(
  'fs_fun' = list(n = number_of_genes)
)


# Assess normalization
# SVM full test sets
svm_full <- assess_normalization(datasets = datasets, 
                                 normalization = normalization, 
                                 fs_fun = fs_fun, 
                                 fit_fun = fit_svm, 
                                 fit_params = svm_params,
                                 pred_fun = predict_svm, 
                                 test_size = NULL, 
                                 n_samplings = 10, 
                                 n_ref_samples = 3, 
                                 ref_class = 'female', 
                                 params = function_params_svm, 
                                 verbose = TRUE, 
                                 requires = 'e1071', 
                                 seed = 1599)

# SVM, max 10 samples in the test set
svm_small <- assess_normalization(datasets = datasets, 
                                  normalization = normalization, 
                                  fs_fun = fs_fun, 
                                  fit_fun = fit_svm, 
                                  fit_params = svm_params,
                                  pred_fun = predict_svm, 
                                  test_size = 10, 
                                  n_samplings = 10, 
                                  n_ref_samples = 3, 
                                  ref_class = 'female', 
                                  params = function_params_svm, 
                                  verbose = TRUE, 
                                  requires = 'e1071', 
                                  seed = 1599)

# Random forest full test sets
rf_full <- assess_normalization(datasets = datasets, 
                                normalization = normalization, 
                                fs_fun = fs_fun, 
                                fit_fun = randomForest, 
                                fit_params = rf_params,
                                pred_fun = predict_rf, 
                                test_size = NULL, 
                                n_samplings = 10, 
                                n_ref_samples = 3, 
                                ref_class = 'female', 
                                params = function_params_rf, 
                                verbose = TRUE, 
                                requires = 'randomForest', 
                                seed = 1599)

# Random forest, max 10 samples in the test set
rf_small <- assess_normalization(datasets = datasets, 
                                 normalization = normalization, 
                                 fs_fun = fs_fun, 
                                 fit_fun = randomForest, 
                                 fit_params = rf_params,
                                 pred_fun = predict_rf, 
                                 test_size = 10, 
                                 n_samplings = 10, 
                                 n_ref_samples = 3, 
                                 ref_class = 'female', 
                                 params = function_params_rf, 
                                 verbose = TRUE, 
                                 requires = 'randomForest', 
                                 seed = 1599)

# kNN full test sets
knn_full <- assess_normalization(datasets = datasets, 
                                 normalization = normalization, 
                                 fs_fun = fs_fun, 
                                 fit_fun = fit_knn, 
                                 fit_params = knn_params,
                                 pred_fun = predict_knn, 
                                 test_size = NULL, 
                                 n_samplings = 10, 
                                 n_ref_samples = 3, 
                                 ref_class = 'female', 
                                 params = function_params_knn, 
                                 verbose = TRUE, 
                                 requires = 'class', 
                                 seed = 1599)

# kNN, max 10 samples in the test set
knn_small <- assess_normalization(datasets = datasets, 
                                  normalization = normalization, 
                                  fs_fun = fs_fun, 
                                  fit_fun = fit_knn, 
                                  fit_params = knn_params,
                                  pred_fun = predict_knn, 
                                  test_size = 10, 
                                  n_samplings = 10, 
                                  n_ref_samples = 3, 
                                  ref_class = 'female', 
                                  params = function_params_knn, 
                                  verbose = TRUE, 
                                  requires = 'class', 
                                  seed = 1599)
saveRDS(object = svm_full, 
        file = 'results/assess_normalization/ratio_a/svm_full.RDS')
saveRDS(object = svm_small, 
        file = 'results/assess_normalization/ratio_a/svm_small.RDS')
saveRDS(object = rf_full, 
        file = 'results/assess_normalization/ratio_a/rf_full.RDS')
saveRDS(object = rf_small, 
        file = 'results/assess_normalization/ratio_a/rf_small.RDS')
saveRDS(object = knn_full, 
        file = 'results/assess_normalization/ratio_a/knn_full.RDS')
saveRDS(object = knn_small, 
        file = 'results/assess_normalization/ratio_a/knn_small.RDS')