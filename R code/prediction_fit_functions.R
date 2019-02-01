# Robin Gradin
#
# Functions for model creation and prediction.
#
# The script contains the following functions:
#   fit_svm
#   predict_svm
#   predict_rf
#   fit_knn
#   predict_knn

library(caret)
library(class)
library(e1071)
library(purrr)
source('R code/helpers.R')

fit_svm <- function(x, y, kernel = 'linear', cost = 1, ...){
  # Fits the SVM model with the linear kernel as default
  # Args:
  #   x: matrix, training data (samples in rows).
  #   y: factor, classes of training data.
  #   kernel; string, SVM kernel.
  #   cost; numeric, SVM cost parameter.
  #   ...: passed to predict.
  # Returns:
  #   Fitted svm model.
  
  model <- svm(x = x, y = y, kernel = kernel, cost = cost, ...)
  return(model)
}

predict_svm <- function(fit, x_test, decision.values = TRUE){
  # Predicts samples in test set using supplied SVM
  # Args:
  #   fit: Fit from fit_svm.
  #   x_test: matrix, samples to be predicted (samples in rows).
  #   decision.values: logical, should decision values be returned.
  # Returns:
  #   Factor with predictions.
  
  pred <- predict(object = fit, 
                  newdata = x_test, 
                  decision.values = decision.values)
  return(pred)
}

predict_rf <- function(fit, x_test, type = 'response'){
  # Predicts test samples using supplied RF model
  # Args:
  #   fit: RandomForest model.
  #   x_test: matrix, samples to be predicted (samples in rows).
  #   type: character, what prediction should be returned?
  # Returns:
  #   Predictions on x_test as factors.
  
  pred <- predict(object = fit, 
                  newdata = x_test, 
                  type = type)
  return(pred)
}

fit_knn <- function(x, y, k = 3){
  # Fit function for knn compatible with the BARA scripts
  # Args:
  #   x: matrix, training data, samples as rows.
  #   y: factor, classes of training data.
  #   k: integer, number of nearest neighbors used for classification.
  # Returns:
  #   List that can be used in prediction functions
  return(list(x = x, y = y, k = k))
}

predict_knn <- function(fit, x_test, prob = FALSE){
  # Classifies samples in x_test using kNN
  # Args:
  #   fit: list, containing model spec, created using fit_knn.
  #   x_test: matrix, samples to be predicted (samples in rows).
  #   prob: logical, should vote probabilities be returned?
  # Returns:
  #   Predictions on x_test as factors.
  
  pred <- knn(
    train = fit$x,
    test = x_test,
    cl = fit$y,
    k = fit$k,
    prob = prob
  )
  return(pred)
}
