# Robin Gradin
#
# Function for tuning model parameters
#
# The script contains the following functions:
#   optimize_parameters

library(doParallel)
library(doRNG)


optimize_parameters <- function(x, 
                                y, 
                                fit_fun, 
                                pred_fun,
                                perf_fun, 
                                model_parameters, 
                                objective = c('maximize', 'minimize'),
                                folds = 10,
                                repeats = 3,
                                keep_predictions = FALSE,
                                cores = 1,
                                seed = 1658,
                                verbose = TRUE,
                                .export = NULL,
                                .packages = NULL){
  # Performs repeated cross-validation over parameter grid.
  #
  # Args:
  #   x: matrix, dataset used for modelling.
  #   y: factor, class labels.
  #   fit_fun: function, fits prediction model.
  #   pred_Fun: function, predicts test samples.
  #   perf_fun: function, calculates performance.
  #   model_parameters: named list, should contain parameter values
  #   objective: character, should the perf_fun be minimized or maximized?
  #   folds: integer, the number of folds used in CV.
  #   repeats: integer, the number of repeats in CV.
  #   keep_predictoins: logical, should the CV-predictions be kept?
  #   cores: integer, number of cores used for calculations.
  #   seed: integer, makes the calculations reproducible.
  #   verbose: logical, should runtime messages be printed?
  #   .export: character, names of objects passed to foreach.
  #   .packages character names of packages passed to foreach
  
  # Initialize variables
  function_call <- match.call()
  objective <- match.arg(arg = objective[1], 
                         choices = c('maximize', 'minimize'))
  set.seed(seed = seed)
  idx_folds <- map(1:repeats, ~createFolds(y = y, k = folds))
  tune_grid <- expand.grid(model_parameters)
  
  # Set up computational environment
  cl <- makeCluster(cores)
  if (cores > 1){
    registerDoParallel(cl)
  } else{
    registerDoSEQ()
  }
  on.exit(expr = {stopCluster(cl)})
  
  # Start the calculations
  results <- foreach (i = 1:nrow(tune_grid), 
                      .export = .export, 
                      .packages = c('purrr', .packages)) %dorng% ({
    # Create lists to store CV output
    performance <- map(1:repeats, ~numeric(folds))
    attr(x = performance, which = 'parameters') <- as.vector(
      x = tune_grid[i, , drop = FALSE],
      mode = 'list'
    )
    
    if (keep_predictions){
      predictions <- map(1:repeats, ~vector(mode = 'list', length = folds))
      attr(x = predictions, which = 'parameters') <- attr(x = performance,
                                                          which = 'parameters')
    } else{
      predictions <- NULL
    }
    
    # Start the repeated cross-validation
    for (j in seq_along(idx_folds)){
      for (k in seq_along(idx_folds[[j]])){
        test_idx <- idx_folds[[j]][[k]]
        fit <- do.call(
          fit_fun, 
          c(
            list(
              x = x[-test_idx, , drop = FALSE], 
              y = y[-test_idx]
            ), 
            as.list(
              tune_grid[i, , drop = FALSE]
            )
          )
        )
        pred <- do.call(
          pred_fun, 
          list(
            fit = fit,
            x_test = x[test_idx, , drop = FALSE]
          )
        )
        performance[[j]][k] <- do.call(
          perf_fun, 
          list(y = y[test_idx], pred = pred)
        )
        if (keep_predictions){
          predictions[[j]][[k]] <- pred
        }
      }
    }
    list(
      performance = performance,
      predictions = predictions
    )
  })
  
  # Find best values in the tuning grid
  perf <- map_dbl(results, ~mean(map_dbl(.$performance, mean)))
  
  if (objective == 'maximize'){
    best_value <- max(perf)
  } else{
    best_value <- min(perf)
  }
  best_idx <- which(perf == best_value)
  best_params <- tune_grid[best_idx, , drop = FALSE]
  
  return(list(best_params = best_params,
              performance = perf,
              results = results,
              tune_grid = tune_grid,
              call = function_call))
}