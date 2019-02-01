# Robin Gradin
#
# Global parameters for the BARA analysis
#
#

svm_params <- list(cost = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50))

rf_params <- list(
  ntree = c(500, 1000, 1500, 2000, 2500, 3000),
  mtry = c(5, 7, 9, 10, 11, 13, 15, 17)
)

knn_params <- list(k = c(1, 2, 3, 4, 5, 6, 7, 8, 9))