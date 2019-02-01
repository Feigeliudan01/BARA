# Robin Gradin
#
# Helper functions.
#
# The script contains the following functions:
#   get_exprs: Returns transposed expression from dataset
#   get_group: Returns group factor from dataset
#   get_name: Returns name of dataset
#   initialize_nested_list: Returns a nested list by recursion
#   create_result_matrix: Returns an empty matrix

get_exprs <- function(dataset, index){
  # Returns the transposed expression matrix of dataset
  # Args:
  #   dataset: list, created at the start of the analysis.
  #   index: integer, which dataset should be returned?
  # Returns:
  #   Transposed expression matrix for specified dataset.
  
  return(t(dataset[[index]]$exprs))
}
get_group <- function(dataset, index){
  # Returns the class factor for dataset at index.
  # Args:
  #   dataset: list, created at the start of the analysis.
  #   index: integer, which dataset should be returned?
  # Returns:
  #   Factor of classes for dataset index.
  
  return(factor(dataset[[index]]$pheno$group))
}
get_name <- function(dataset, index){
  # Returns the name of dataset at index.
  # Args:
  #   dataset: list, created at the start of the analysis.
  #   index: integer, which dataset should be returned?
  # Returns:
  #   Character of the dataset's name.
  
  data_name <- unique(dataset[[index]]$pheno$dataset)
  if (length(data_name) > 1){
    stop('More than 1 name was found.')
  }
  return(data_name)
}

initialize_nested_list <- function(size){
  # Create an empty nested list.
  # Args: 
  #   size: numeric, specifiying the number of levels (length(size))
  #         and the number of lists in each level. ex c(2, 3, 4).
  # Returns:
  #   Empty nested list.
  
  if (length(size) == 1){
    return(vector(mode = 'list', length = size))
  } else{
    map(1:size[1], ~initialize_nested_list(size[-1]))
  }
}
create_result_matrix <- function(dataset){
  # Create empty matrix for storing prediction results.
  # Args:
  #   dataset: list, created at the start of the analysis.
  # Returns:
  #   Empty matrix of size nxn, where n = length(dataset)
  
  n <- length(dataset)
  sample_names = names(dataset)
  
  result_matrix <- matrix(
    data = NA, 
    nrow = n, 
    ncol = n, 
    dimnames = list(
      'Training set' = sample_names,
      'Test set' = sample_names
      )
  )
  return(result_matrix)
}
get_empty_performance_list <- function(dataset){
  # Create list of result matrices.
  # Args:
  #   dataset: list, created at the start of the analysis.
  # Returns 
  #   Named list with empty result matrices (accuracy, cross_entropy, mcc)
  
  perf <- list(
    'accuracy' = list(
      'lm' = create_result_matrix(dataset)
    ),
    'cross_entropy' = list(
      'lm' = create_result_matrix(dataset)
    ),
    'mcc' = list(
      'lm' = create_result_matrix(dataset)
    )
  )
  return(perf)
}