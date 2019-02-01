# Robin Gradin
#
# Functions for feature selection.
#
# The script contains the following functions:
#   select_significant_limma

library(limma)

select_significant_limma <- function(x, y, n){
  # Fits limma to data and returns n most significant genes.
  # x: matrix, expression data.
  # y: factor, groups or true classes.
  # n: integer, number of genes to return.
  #
  # returns:
  #   selected_genes: character of length n, identified genes.
  #   fit: moderated limma object.
  
  model <- model.matrix(~y)
  fit <- lmFit(object = t(x), design = model)
  fit <- eBayes(fit)
  selected_genes <- rownames(topTable(fit, coef = 2, number = n))
  return(list(
    'selected_genes' = selected_genes,
    'fit' = fit
  ))
}
