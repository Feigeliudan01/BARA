# Robin Gradin
#
# Batch adjustment by reference alignment (2019)
# Prediction results at different number of 
# reference samples for BARA

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

# Normalization fit function
# Estimates compression and 
# returns a compressed training set.
bara_fit_wrapper <- function(x, ref, ...){
  input_vars <- match.call()
  fit <- input_vars$project_train(x = x,
                                n_dimensions = NULL,
                                loss = 0.1, 
                                scale_var = FALSE)
  fit$loss <- 0.1
  fit$x_org <- x
  x <- reconstruct_data(object = fit)
  return(list(
    x_train = x,
    fit = fit,
    ref = ref
  ))
}

normalization <- list(
  fit = bara_fit_wrapper,
  transform = function(fit, x_test, ref){
    nrm <- bara(x_train = fit$fit$x_org,
                x_test = x_test, 
                reference_train = fit$ref, 
                reference_test = ref, 
                n_dimensions = NULL,
                loss = fit$fit$loss, 
                scale_var = FALSE)
    return(list(x_test = nrm$x_test, ref = ref))
  }
)

function_params_svm <- list(
  'fs_fun' = list('n' = number_of_genes),
  'fit_fun' = list('kernel' = 'linear'),
  'pred_fun' = list('decision.values' = FALSE),
  'normalization.fit' = list('fit_fun' = fit_svm,
                             'pred_fun' = predict_svm,
                             'requires' = 'e1071',
                             'project_train' = project_train)
)
function_params_rf <- list(
  'fs_fun' = list('n' = number_of_genes),
  'normalization.fit' = list('fit_fun' = randomForest,
                             'pred_fun' = predict_rf,
                             'requires' = 'randomForest',
                             'project_train' = project_train)
)

function_params_knn <- list(
  'fs_fun' = list('n' = number_of_genes),
  'normalization.fit' = list('fit_fun' = fit_knn,
                             'pred_fun' = predict_knn,
                             'requires' = 'class',
                             'project_train' = project_train)
)


# Assess normalization
n_ref <- 1:12
for (i in seq_along(n_ref)){
  # SVM full test sets
  svm_full <- assess_normalization(datasets = datasets, 
                                   normalization = normalization, 
                                   fs_fun = fs_fun, 
                                   fit_fun = fit_svm, 
                                   fit_params = svm_params,
                                   pred_fun = predict_svm, 
                                   test_size = NULL, 
                                   n_samplings = 10, 
                                   n_ref_samples = i,
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
                                    n_ref_samples = i,
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
                                  n_ref_samples = i,
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
                                   n_ref_samples = i,
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
                                   n_ref_samples = i,
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
                                    n_ref_samples = i,
                                    ref_class = 'female', 
                                    params = function_params_knn, 
                                    verbose = TRUE, 
                                    requires = 'class', 
                                    seed = 1599)
  saveRDS(object = svm_full, 
          file = sprintf('results/bara_ref_samples/svm_full_%s.RDS', i)
          )
  saveRDS(object = svm_small, 
          file = sprintf('results/bara_ref_samples/svm_small_%s.RDS', i)
          )
  saveRDS(object = rf_full, 
          file = sprintf('results/bara_ref_samples/rf_full_%s.RDS', i)
          )
  saveRDS(object = rf_small, 
          file = sprintf('results/bara_ref_samples/rf_small_%s.RDS', i)
          )
  saveRDS(object = knn_full, 
          file = sprintf('results/bara_ref_samples/knn_full_%s.RDS', i)
          )
  saveRDS(object = knn_small, 
          file = sprintf('results/bara_ref_samples/knn_small_%s.RDS', i)
          )
}
