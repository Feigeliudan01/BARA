# Robin Gradin
# 
# Runs the scripts for all normalization methods

# Load the paths to the scripts
scripts <- list.files(path = 'src/assess_normalization/',
                      pattern = '.R$', 
                      full.names = TRUE)
this_script <- which(scripts == 'src/assess_normalization/run_all_scripts.R')
scripts <- scripts[-this_script]

# Run all scripts
for (i in seq_along(scripts)){
  eval_env <- new.env()
  writeLines(sprintf(
    '%s| %s',
    Sys.time(),
    scripts[i]
  ))
  log_file <- sprintf('%s/Log_%s.txt',
                      dirname(scripts[i]),
                      basename(scripts[i])
  )
  output <- file(log_file)
  sink(file = output)
  try(
    expr = {
      sys.source(file = scripts[i],
                 envir = eval_env, 
                 toplevel.env = eval_env)
    }, 
    outFile = log_file
  )
  sink(file = NULL)
  close(output)
  rm(eval_env)
}
