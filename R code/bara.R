# Robin Gradin
#
# The BARA function implemented in manuscript
#
# This script contains the following functions:
#   bara
#   create_logical_vector
#   get_compression

library(caret)
library(doParallel)
library(doRNG)
library(matrixStats)
library(purrr)

bara <- function(x_train,
                 x_test,
                 reference_train,
                 reference_test,
                 batches = NULL,
                 n_dimensions = NULL,
                 loss = 0.01,
                 scale_var = FALSE){
  # Adjusts batch effects between datasets for prediction problems.
  # Args:
  #   x_train: matrix, training data (samples as rows).
  #   x_test: matrix, test data (samples as rows).
  #   reference_train: logical or numeric, describes the indices
  #                    of the reference samples in the training set.
  #   reference_test: logical or numeric, describes the indices
  #                   of the reference samples in the test set.
  #   batches: character, describing batches (if any) in the test set.
  #   n_dimensions: integer, number of dimensions retained during compression.
  #   loss: numeric, can be used instead of n_dimensions to describe the
  #         amount of compression applied.
  #   scale_var: logical, should the variance be standardized?
  # Returns:
  #   x_train: matrix, compressed training set.
  #   x_test: matrix, normalized test set.
  #   n_dimensions: integer, number of dimensions retained during compression.
  #   loss: numeric, fraction of variance in x_train lost during compression.
  
  # Check input data
  if (ncol(x_train) != ncol(x_test)){
    stop(
      sprintf(
        'x_train (%s) and x_test (%s) should have the same number of columns.',
        ncol(x_train),
        ncol(x_test)
      )
    )
  }
  if (any(colnames(x_train) != colnames(x_test))){
    stop('The column names of x_train and x_test must be match')
  }
  if (is.null(batches)){
    batches <- rep('one_batch', nrow(x_test))
  }
  if(!all(unique(batches) %in% batches[reference_test])){
    stop('At least one reference sample must be present in each batch.')
  }
  if(!is.matrix(x_train)){
    x_train <- as.matrix(x_train)
    message('x_train was cast to matrix.')
  }
  if(!is.matrix(x_test)){
    x_test <- as.matrix(x_test)
    message('x_test was cast to matrix.')
  }
  if (is.numeric(reference_train)){
    reference_train <- create_logical_vector(
      size = nrow(x_train),
      idx = reference_train
    )
  }
  if (is.numeric(reference_test)){
    reference_test <- create_logical_vector(
      size = nrow(x_test),
      idx = reference_test
    )
  }
  
  # Save varnames for later
  # and get unique batches
  varnames <- colnames(x_train)
  unique_batches <- unique(batches)
  
  # Subtract the mean vector from the datasets
  mean_vector <- colMeans(x_train)
  x_train <- sweep(x = x_train, MARGIN = 2, STATS = mean_vector, FUN = '-')
  x_test <- sweep(x = x_test, MARGIN = 2, STATS = mean_vector, FUN = '-')
  
  # Scale if desired
  if (scale_var){
    sd_vector <- colSds(x_train)
    x_train <- sweep(x = x_train, MARGIN = 2, STATS = sd_vector, FUN = '/')
    x_test <- sweep(x = x_test, MARGIN = 2, STATS = sd_vector, FUN = '/')
  }
  
  # Decompose the training data
  usv <- svd(x_train)
  captured_variance <- cumsum((usv$d^2) / sum(usv$d^2))
  
  # Determine compression
  dimensions <- get_compression(
    captured_variance = captured_variance,
    n_dimensions = n_dimensions,
    loss = loss,
    min_dim = min(dim(x_train))
  )
  compression_factor <- captured_variance[dimensions]
  
  # Project onto right singular vector
  prj_train <- x_train %*% usv$v[, 1:dimensions, drop = FALSE]
  prj_test <- x_test %*% usv$v[, 1:dimensions, drop = FALSE]
  
  # Get position of train's reference samples
  reference_position <- colMeans(prj_train[reference_train, , drop = FALSE])
  
  # Adjust all batches
  for (i in seq_along(unique_batches)){
    idx_batch <- batches == unique_batches[i]
    test_position <- colMeans(
      prj_test[idx_batch & reference_test, , drop = FALSE]
      )
    correction <- reference_position - test_position
    prj_test[idx_batch, ] <- sweep(
      x = prj_test[idx_batch, , drop = FALSE],
      MARGIN = 2,
      STATS = correction,
      FUN = '+'
    )
  }
  
  # Reconstruct the data
  x_train <- prj_train %*% t(usv$v[, 1:dimensions, drop = FALSE])
  x_test <- prj_test %*% t(usv$v[, 1:dimensions, drop = FALSE])
  
  # Scale
  if (scale_var){
    x_train <- sweep(x = x_train, MARGIN = 2, STATS = sd_vector, FUN = '*')
    x_test <- sweep(x = x_test, MARGIN = 2, STATS = sd_vector, FUN = '*')
  }
  
  # Add mean expression
  x_train <- sweep(x = x_train, MARGIN = 2, STATS = mean_vector, FUN = '+')
  x_test <- sweep(x = x_test, MARGIN = 2, STATS = mean_vector, FUN = '+')
  
  # Add variable names
  colnames(x_train) <- varnames
  colnames(x_test) <- varnames
  
  bara_results <- list(
    x_train = x_train,
    x_test = x_test,
    n_dimensions = dimensions,
    loss = 1 - compression_factor
  )
  return(bara_results)
}

create_logical_vector <- function(size, idx){
  # Creates a logical vector with TRUE at idx.
  # Args: 
  #   size: numeric, length of the final vector
  #   idx: numeric, which indices should contain TRUE.
  # Returns:
  #   Logical vector with TRUE only at indices specified by idx.
  
  log_vec <- rep(FALSE, size)
  log_vec[idx] <- TRUE
  return(log_vec)
}

get_compression <- function(captured_variance, 
                            n_dimensions, 
                            loss,
                            min_dim){
  # Determines the amount of compression to apply in BARA.
  # Args:
  #   captured_variance: numeric, cumulated sum of amount of variance
  #                      captured by each singular vector.
  #   n_dimensions: numeric, number of dimensions to be retained.
  #   loss: numeric, amount of variance that are acceptable to lose.
  #   min_dim, minimum number of dimensions possible to retain.
  # Returns:
  #   Integer describing the number of dimensions to retain.
  if (is.null(n_dimensions)){
    if (loss <= 0 || loss > 1){
      stop('Loss must be between 0 and 1')
    } else{
      dimensions <- which((1 - captured_variance) < loss)[1]
    }
  } else{
    if (n_dimensions > min_dim || n_dimensions < 1){
      stop(
        sprintf(
          'n_dimensions must be between 1 and %s.',
          min_dim
        )
      )
    } else{
      dimensions <- n_dimensions
    }
  }
  return(dimensions)
}

cv_bara <- function(x, y,
                    fit_fun,
                    pred_fun,
                    perf_fun,
                    objective = c('maximize', 'minimize'),
                    folds = 10,
                    repeats = 3,
                    frac = 0,
                    scale_var = FALSE,
                    cores = 1,
                    .export = NULL,
                    .packages = NULL,
                    seed = 82347,
                    verbose = FALSE){
  # Cross-validation (CV) to determine compression in BARA
  # Args:
  #   x: matrix, training data (samples as rows).
  #   y: factor, true classes of training data.
  #   fit_fun: function, takes x and y as input and returns a prediction model.
  #   pred_fun: function, takes fit and x_test as input and 
  #             generates predictions on x_test.
  #   perf_fun: function, takes y (true class) and pred as input, 
  #             and returns a numeric performance metric. 
  #   objective: character, should the performance metric 
  #              be maximized or minimized?
  #   folds: integer, number of folds in CV.
  #   repeats: integer, number of repeats in repeated CV.
  #   frac: numeric, fraction of dimensions to remove in each iteration.
  #         All dimensions are examined if frac is set to 0.
  #   scale_var: logical, should the variance be standardized?
  #   cores: integer, number of cores to use in parallel computations.
  #   .export: character, passed to foreach.
  #   .packages: character, passed to foreach.
  #   seed: integer, to make results reproducible.
  #   verbose: logical, should runtime messages be displayed.
  # Returns:
  #

  # Check input
  if (nrow(x) != length(y)){
    stop('The number of rows in x must match the length of y')
  }
  objective <- match.arg(
    arg = objective[1],
    choices = c('maximize', 'minimize')
  )
  if (frac < 0 || frac >=1){
    stop('frac must be a number between 0 and 1.')
  }
  if (!is.matrix(x)){
    x <- as.matrix(x)
    message('x was cast to matrix.')  # Not affected by verbose
  }
  # Remove zero-variance columns
  zero_variance <- which(colVars(x) == 0)
  if (length(zero_variance) > 1){
    x <- x[, -zero_variance, drop = FALSE]
    message(
      sprintf(
        'Zero-variance columns were identified, %s column/s were removed',
        length(zero_variance)
      )
    )
  }

  # Get the CV-splits
  set.seed(seed)
  idx_folds <- map(1:repeats, ~createFolds(y = y, k = folds))

  # Calculate the maximum number of dimensions evaluated
  max_dimension <- min(
    (nrow(x) - ceiling(nrow(x) / folds)),
    ncol(x)
  )

  # Calculate which dimensions should be examined
  dim_examined <- get_dimensions(max_dimension = max_dimension, frac = frac)
  if (verbose){
    message(sprintf('Analyzing %s dimensions', length(dim_examined)))
  }

  # Setup computational environment
  cl <- makeCluster(cores)
  if (cores > 1){
    registerDoParallel(cl)
    if (verbose){
      message('Runtime messages are suppressed during parallel computations.')
    }
  } else{
    registerDoSEQ()
  }
  on.exit(expr = {stopCluster(cl)})

  results <- foreach(i = seq_along(dim_examined),
                     .export = c('bara_cv_compress', .export),
                     .packages = c('purrr', .packages)) %dorng% ({
      performance <- map(1:repeats, ~numeric(folds))
      attr(performance, which = 'dimensions') <- dim_examined[i]

      if (verbose && cores == 1){
        message(
          sprintf(
            '%s| Retaining %s singular vector/s.',
            Sys.time(),
            dim_examined[i]
          )
        )
      }

      # Start repeated CV
      for (j in 1:repeats){
        for (k in 1:folds){
          test_idx <- idx_folds[[j]][[k]]
          # Run BARA compression
          bara_data <- bara_cv_compress(
            x_train = x[-test_idx, , drop = FALSE],
            x_test = x[test_idx, , drop = FALSE],
            dimensions = dim_examined[i],
            scale_var = scale_var
          )
          # Build model
          fit <- do.call(
            fit_fun,
            list(
                x = bara_data$x_train,
                y = y[-test_idx]
              )
          )
          # Classify test set
          pred <- do.call(
            pred_fun,
            list(
              fit = fit,
              x_test = bara_data$x_test
            )
          )
          # Evaluate performance
          perf <- do.call(
            perf_fun,
            list(
              y = y[test_idx],
              pred = pred
            )
          )
          performance[[j]][[k]] <- perf
        }
      }
      performance
    }
  )

  # Summarize the results
  perf <- map_dbl(results, ~mean(map_dbl(., mean)))
  attr(x = perf, which = 'dimensions') <- map_dbl(results, attr, 'dimensions')

  if (objective == 'maximize'){
    best_value <- max(perf)
  } else{
    best_value <- min(perf)
  }
  best_idx <- which(perf == best_value)
  best_dim <- min(attr(x = perf, which = 'dimensions')[best_idx])

  cv_bara_results <- list(
    best_dimension = best_dim,
    performance = perf,
    cv_results = results
  )
  
  return(cv_bara_results)
}


get_dimensions <- function(max_dimension, frac){  
  # Get the dimensions evaulated in BARA CV
  # Args:
  #   max_dimensions: integer, maximum number of dimensions to be evaluated.
  #   frac: fraction of dimensions to be skipped in each iteration.
  # Returns:
  #   Integer vector, the dimensions that should be evaluated.
  if (frac == 0){
    dimensions <- seq(from = 1, to = max_dimension)    
  } else{
    dimensions <- c(max_dimension)
    while(dimensions[1] > 1){
      n_to_rm <- floor(dimensions[1] * frac)
      if (n_to_rm < 1){
        n_to_rm <- 1
      }           
      dimensions <- c(
        dimensions[1] - n_to_rm,
        dimensions
      )
    }
  }
  return(dimensions)
}

bara_cv_compress <- function(x_train, x_test, dimensions, scale_var){
  # Perform the BARA compression without adjusting for batch effects.
  # Args: 
  #   x_train: matrix, internal train set.
  #   x_test: matrix, internal test set.
  #   dimensions: integer, number of dimensions to retain.
  #   scale_var: logical, should the variance be standardized?
  # Returns:
  #   x_train: matrix, compressed training set.
  #   x_test: matrix, test set projected onto training sets space.

  var_names <- colnames(x_train)
  mean_vector <- colMeans(x_train)
  x_train <- sweep(x = x_train, MARGIN = 2, STATS = mean_vector, FUN = '-')
  x_test <- sweep(x = x_test, MARGIN = 2, STATS = mean_vector, FUN = '-')

  if (scale_var){
    sd_vector <- colSds(x_train)
    x_train <- sweep(x = x_train, MARGIN = 2, STATS = sd_vector, FUN = '/')
    x_test <- sweep(x = x_test, MARGIN = 2, STATS = sd_vector, FUN = '/')
  }
  
  # Compress and reconstruct the datasets
  usv <- svd(x_train, nv = dimensions)
  if (ncol(usv$v) != dimensions){
    stop(
      'Something went wrong with svd. Unexpected number of singular vectors.'
    )
  }
  x_train <- (x_train %*% usv$v) %*% t(usv$v)
  x_test <- (x_test %*% usv$v) %*% t(usv$v)

  if (scale_var){
    x_train <- sweep(x = x_train, MARGIN = 2, STATS = sd_vector, FUN = '*')
    x_test <- sweep(x = x_test, MARGIN = 2, STATS = sd_vector, FUN = '*')
  }

  x_train <- sweep(x = x_train, MARGIN = 2, STATS = mean_vector, FUN = '+')
  x_test <- sweep(x = x_test, MARGIN = 2, STATS = mean_vector, FUN = '+')
  colnames(x_train) <- var_names
  colnames(x_test) <- var_names
  
  return(list(
    x_train = x_train,
    x_test = x_test
  ))
}

project_train <- function(x, 
                           n_dimensions = NULL,
                           loss = 0.01,
                           scale_var = FALSE){
  # Calculates the projection of x onto its singular vector.
  # Args:
  #   x: matrix, dataset to be projected onto its singular vectors.
  #   n_dimensions: numeric, number of dimensions to retain.
  #   loss: numeric, acceptable fraction of variance to lose 
  #           during compression. Evaluated if n_dimensions = NULL.
  #   scale_var: logical, should the vairance be standardized?
  # Returns:
  #   x: matrix, projection onto the singular vectors.
  #   singular_vectors: matrix, the right singular vectors.
  #   mean_vector: numeric, original column means,
  #   sd_vector: numeric, original column standard deviations.
  #   var_names: character, column names of the input matrix.
  var_names <- colnames(x)
  x <- scale(x, center = TRUE, scale = scale_var)
  mean_vector <- attr(x = x, which = 'scaled:center')
  if(scale_var){
    sd_vector <- attr(x = x, which = 'scaled:scale')
  } else{
    sd_vector = FALSE
  }
  # Decompose the data
  usv <- svd(x)
  captured_variance <- cumsum((usv$d^2) / sum(usv$d^2))
  n_dim <- get_compression(
    captured_variance = captured_variance,
    n_dimensions = n_dimensions,
    loss = loss,
    min_dim = min(dim(x))
  )

  x <- x %*% usv$v[, 1:n_dim, drop = FALSE]  

  return(list(
    x = x,
    singular_vectors = usv$v[, 1:n_dim, drop = FALSE],
    mean_vector = mean_vector,
    sd_vector = sd_vector,
    var_names = var_names
  ))
}

project_test <- function(object, x){
  # Calculates the projection of x onto the training sets singular vectors.
  # Args:
  #   object: list, output from project_train().
  #   x: matrix, test set to be projected.
  # Returns:
  #   Matrix, projection onto the singular vectors

  # Check input
  if (!all(colnames(x) == object$var_names)){
    stop('x must have the same column names as the training set.')
  }
  # Standardize
  x <- scale(x, center = object$mean_vector, scale = object$sd_vector)
  # Project
  x <- x %*% object$singular_vectors
  return(x)
}

reconstruct_data <- function(object, x = NULL){
  # Reconstruct the data to the original dimensions.
  # Args:
  #   object: list, output from project_train().
  #   x: optional matrix, data to be reconstructed. If NULL, 
  #      the data is taken from object.
  # Returns:
  #   Matrix, reconstructed data.

  if (is.null(x)){
    x <- object$x %*% t(object$singular_vectors)
  } else{
    x <- x %*% t(object$singular_vectors)
  }
  if (!is.logical(object$sd_vector)){
    x <- sweep(x = x, MARGIN = 2, STATS = object$sd_vector, FUN = '*')
  }
  x <- sweep(x = x, MARGIN = 2, STATS = object$mean_vector, FUN = '+')
  colnames(x) <- object$var_names
  return(x)
}