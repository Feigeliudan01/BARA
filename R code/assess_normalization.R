# Robin Gradin
#
# Batch adjustment by reference alignment (2019)
# Function to assess the performance after normalization

# Load dependencies
library(bapred)
library(ggplot2)
library(ggthemes)
library(randomForest)
library(reshape2)
library(RGvisualization)
library(stringr)
library(viridis)

# In this script
call_function <- function(FUN, input, params, name = 'fs_fun'){
  if (name %in% names(params)){
    input <- c(input, params[[name]])
  }
  function_output <- do.call(
    what = FUN,
    args = input
  )
  return(function_output)
}


assess_normalization <- function(datasets,
                                 normalization,
                                 fs_fun,
                                 fit_fun,
                                 fit_params,
                                 pred_fun,
                                 test_size = NULL,
                                 n_samplings = 1,
                                 n_ref_samples = 3,
                                 ref_class = NULL,
                                 params = list(),
                                 verbose = TRUE,
                                 requires = c(),
                                 seed = 3276){
  # Runs normalization and calculates performance metrics.
  # Args:
  #   datasets:      List, contains the datasets.
  #   normalization: Named list of functions. Should contain fit and transform.
  #                  fit is applied to the training set and transform
  #                  to the test set. fit must take x, ref and ... as input, 
  #                  and transform must take fit, x_test and ref as input. 
  #                  fit should return a list with x_train, fit and ref.
  #                  transform should return a list with x_test.
  #   fs_fun:        Function, takes x (matrix) and y (factor) as input
  #                  and returns a character vector with selected genes.
  #   fit_fun:       Function, must take x (matrix training set) and 
  #                  y (factor, class training set) as input.
  #   fit_params:    Named list of parameters used to tune the fit function.
  #   pred_fun:      Function, must take fit (output from fit_fun),
  #                  and x_test (matrix, test set) as input.
  #   test_size:     Numeric or NULL, desired size of the test set.
  #                  If NULL or test_size > test set, all samples in the
  #                  test set will be used for predictions.
  #   n_samplings:   Integer, describes the number of times that the test
  #                  set is samples for test size. For example, 
  #                  if test_size=6 and n_samples=5, 6 samples are randomly
  #                  selected from the test set and normalized and predicted.
  #                  This is repeated 5 times.
  #   n_ref_samples: Integer, describes the number of samples 
  #                  reserved as refence samples. These are included in the
  #                  test size. I.e. no more samples than test_size are used
  #                  as test set. The performance is not calculated on these
  #                  samples.
  #   ref_class:     Character, describes the class that should be sampled
  #                  for reference samples.
  #   params:        List of parameters sent to the different functions.
  #                  For example, list(fit_fun = list(cost = 1)) passes a
  #                  value of 1 for the cost parameter in the fit fun.
  #   verbose:       Logical, should runtime messages be displayed?
  #   requires:      Character, required packages to run all functions.
  # Returns:
  #   Nested list. First level describes the training set, second level
  #   contains the test sets, and the third level contains the samplings
  #   on the test set.
  set.seed(seed = seed)
  
  if (!is.null(test_size)){
    if (n_ref_samples >= test_size){
      stop(
        sprintf(
          'n_ref_samples (%s) should be smaller than test_size (%s).',
          n_ref_samples,
          test_size
        )
      )
    }
  }
  n_ref_samples_input <- n_ref_samples
  
  if (verbose){
    writeLines(
      sprintf(
        '%s| Starting assess_normalization().',
        Sys.time()
      )
    )
  }
  # Sequences to iterate through
  iteration = list(first = seq_along(datasets),
                   third = seq(from = 1, to = n_samplings))
  
  # Just select the first level in the first dataset if empty
  if (is.null(ref_class)){
    ref_class = levels(get_group(datasets, 1))[1]
  }
  
  # Start the iterations
  result_list <- list()
  for (i in iteration$first){
    
    if (verbose){
      writeLines(
        sprintf(
          '%s| Dataset %s of %s.',
          Sys.time(),
          i,
          length(iteration$first)
        )
      )
    }
    
    # Get training set
    x_train <- get_exprs(datasets, i)
    y_train <- get_group(datasets, i)
    train_name <- get_name(datasets, i)
    result_list[[train_name]] <- list()
    
    # Select reference samples from the training set
    ref_train <- sample(x = which(y_train == ref_class), 
                        size = n_ref_samples, 
                        replace = FALSE)
    
    # Find gene signature
    genes <- call_function(
      FUN = fs_fun,
      input = list(x = x_train, y = y_train),
      param = params,
      name = 'fs_fun'
    )
    
    # Normalize training set
    norm_fit <- call_function(
      FUN = normalization$fit,
      input = list(x = x_train[, genes, drop = FALSE],
                   ref = ref_train,
                   y = y_train),
      params = params,
      name = 'normalization.fit'
    )
    
    # Tune prediction model
    fit_tune <- optimize_parameters(x = norm_fit$x_train, 
                                    y = y_train, 
                                    fit_fun = fit_fun,
                                    pred_fun = pred_fun, 
                                    perf_fun = performance_mcc, 
                                    model_parameters = fit_params,
                                    objective = 'maximize', 
                                    folds = 10, 
                                    repeats = 3, 
                                    keep_predictions = FALSE, 
                                    cores = 8, 
                                    seed = 1682, 
                                    verbose = FALSE,
                                    .export = c(),
                                    .packages = requires
                                    )
    # Define prediction model
    fit_model <- call_function(
      FUN = fit_fun,
      input = c(
                list(
                     x = norm_fit$x_train,
                     y = y_train
                     ),
                as.vector(
                     x = fit_tune$best_params[1, , drop = FALSE],
                     mode = 'list'
                     )
                     ), 
      params = params, 
      name = 'fit_fun'
    )
    
    # Iterate over the test set
    test_set_iteration <- iteration$first[-i]
    for (j in test_set_iteration){
      # Extract test data
      x_test <- get_exprs(datasets, j)
      y_test <- get_group(datasets, j)
      test_name <- get_name(datasets, j)
      result_list[[train_name]][[test_name]] <- list()
      
      # Make sure that the test set sizes are ok.
      # Multiple sampling is still ok to assess impact of normalization
      if(is.null(test_size)){
        test_size <- nrow(x_test)
      } else if (test_size > nrow(x_test)){
        test_size = nrow(x_test)
      }
      
      if (n_ref_samples >= test_size){
        warning(
          sprintf(
            'n_ref_samples was decreased for dataset %s',
            j
          )
        )
        n_ref_samples <- test_size - 1
      }
      
      # Start the inner prediction loop
      for (k in iteration$third){
        # Select reference samples
        ref_test <- sample(x = which(y_test == ref_class),
                           size = n_ref_samples,
                           replace = FALSE)
        selected_tests <- sample(x = seq(from = 1, 
                                         to = nrow(x_test)-length(ref_test), 
                                         by = 1),
                                 size = (test_size - n_ref_samples),
                                 replace = FALSE)
        selected_tests <- c(ref_test, selected_tests)
        ref_test <- seq(from = 1, to = n_ref_samples, by = 1)
        
        
        # Normalize test set
        norm_transform <- call_function(
          FUN = normalization$transform,
          input = list(fit = norm_fit, 
                       x_test = x_test[selected_tests, genes, drop = FALSE],
                       ref = ref_test),
          params = params,
          name = 'normalization.transform'
        )
        
        # Predict samples
        pred <- call_function(
          FUN = pred_fun,
          input = list(
            fit = fit_model,
            x_test = norm_transform$x_test[-ref_test, , drop = FALSE]
          ),
          params = params,
          name = 'pred_fun'
        )
        
        # Calculate performance
        perf <- performance_mcc(y = y_test[selected_tests][-ref_test],
                                pred = pred,
                                pos = y_train[1])
        result_list[[train_name]][[test_name]][[k]] <- perf
        
        # Restore this value if it was altered
        n_ref_samples <- n_ref_samples_input
      }
    }
  }
  return(result_list)
}
