# Robin Gradin
#
# Functions for calculating performance metrics.
#
# The script contains the following functions:
#   performance_accuracy
#   performance_svm_CE
#   performance_mcc


performance_accuracy = function(y, pred){
  # Calculates prediction accuracy.
  # Args:
  #   y: factor, true classes.
  #   pred: factor, predicted classes.
  # Returns:
  #   Prediction accuracy.
  
  accuracy = mean(y == pred) * 100
  return(accuracy)
}

performance_svm_CE <- function(y, pred){
  # Calculates the mean cross-entropy error of predictions.
  # Args:
  #   y: factor, true classes.
  #   pred: factor, predicted classes with decision values
  #         as attribute.
  # Returns:
  #   Mean cross-entropy error
  
  dval <- abs(as.numeric(attr(x = pred, which = 'decision.values')))
  misclassified <- y != pred
  dval[misclassified] <- dval[misclassified] * -1
  sigmoid_dval <- 1 / (1 + exp(-dval))
  mean_cross_entropy <- -mean(log(sigmoid_dval))
  return(mean_cross_entropy)
}

performance_mcc <- function(y, pred, pos = levels(y)[2]){
  # Calculates Mathews correlation coefficient of predictions.
  # Args:
  #   y: factor, true classes.
  #   pred: factor, predicted classes.
  # Returns:
  #   Mathews correlation coefficient.
  
  if (nlevels(y) > 2){
    stop(sprintf('The response variable should have 2 levels, %s supplied.', 
                 nlevels(y)))
  }
  tp <- sum((y == pred) & (y == pos))
  tn <- sum((y == pred) & (y != pos))
  fp <- sum((y != pred) & (pred == pos))
  fn <- sum((y != pred) & (pred != pos))
  denominator <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
  if (denominator == 0){
    denominator <- 1
  }
  mcc <- ((tp * tn) - (fp * fn)) / sqrt(
    denominator
  )
  return(mcc)
}
