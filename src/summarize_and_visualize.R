# Robin Gradin
#
# Summarize and visualize the normalization results
#

library(data.table)
library(ggplot2)
library(magrittr)
library(purrr)
library(reshape2)
library(stringr)

name_changer <- c(
  'bara' = 'BARA',
  'mean_centering' = 'Mean centered',
  'ratio_a' = 'Ratio A',
  'none' = 'None',
  'combat' = 'Combat',
  'fabatch' = 'FAbatch',
  'fsva_fast' = 'fSVA fast',
  'fsva_exact' = 'fSVA exact',
  'ref_centered' = 'ref-Centered',
  'ref_ratio' = 'ref-Ratio A',
  'standardized' = 'Standardized'
)

result_files <- list.files(path = 'results/assess_normalization/',
                           pattern = '.RDS$', 
                           recursive = TRUE, 
                           full.names = TRUE)
results <- list()
writeLines(sprintf('%s result files', length(result_files)))
for (i in seq_along(result_files)){
  path_parts <- str_split(string = result_files[[i]],
                          pattern = '/')
  norm_method <- path_parts[[1]][4]
  model <- path_parts[[1]][5] %>% gsub(pattern = '.RDS$', 
                                       replacement = '')
  results[[model]][[name_changer[norm_method]]] <- readRDS(
    file = result_files[i]
  )
}

# Summarize the results for each test set
sum_test_set <- map(results, function(mdl){
  map(mdl, function(norm_method){
    map(norm_method, function(training_set){
      map(training_set, ~mean(unlist(.)))
    })
  })
})

# Summarize the results for each training set
sum_training_set <- map(sum_test_set, function(mdl){
  map(mdl, function(norm_method){
    map_dbl(norm_method, ~median(unlist(.)))
  })
})

# Create matrices for each model
result_matrices <- map(sum_training_set, function(mdl){
  as.data.frame(t(do.call(
    what = rbind,
    args = mdl
  )))
})

mean_mcc <- map(result_matrices, function(res){
  data.frame(
    'Mean' = round(apply(res, 2, mean), digits = 3),
    'Stdev' = round(apply(res, 2, sd), digits = 3)
  )
})
for (i in seq_along(mean_mcc)){
  mean_mcc[[i]]$output <- sprintf('%.2f \u00B1 %.2f', mean_mcc[[i]]$Mean, mean_mcc[[i]]$Stdev)
}
result_matrices <- map(result_matrices, melt, id.vars = NULL)

# Define plot function
result_boxplot <- function(x, title = NULL, subtitle = NULL){
  gg <- ggplot(data = x, aes(x = variable, y = value)) +
  geom_boxplot(aes(colour = variable, fill = variable), width = 0.75) +
  stat_summary(
    fun.y = 'median',
    mapping = aes(ymax = ..y.., ymin = ..y..),
    geom = 'crossbar',
    width = 0.5,
    fatten = 0,
    colour = 'white'
  ) +
  scale_color_grey() +
  scale_fill_grey() +
  theme(
    panel.background = ggplot2::element_rect(fill = NA),
    legend.background = ggplot2::element_rect(fill = NA, colour = NA),
    legend.box.background = ggplot2::element_rect(fill = NA, colour = NA),
    plot.subtitle = element_text(colour = 'dimgrey'),
    legend.key = ggplot2::element_rect(fill = NA, colour = NA),
    axis.line = ggplot2::element_line(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  ylim(c(-1, 1)) +
  labs(
    x = NULL,
    y = 'Matthews Correlation Coefficient',
    colour = 'Normalization',
    fill = 'Normalization',
    title = title,
    subtitle = subtitle
  )
  return(gg)
}

box_plots <- list()
box_plots[['knn_full']] <- result_boxplot(x = result_matrices$knn_full,
               title = 'Prediction Performance',
               subtitle = 'kNN')

box_plots[['knn_small']] <- result_boxplot(x = result_matrices$knn_small,
               title = 'Prediction Performance',
               subtitle = 'kNN, small test sets')

box_plots[['rf_full']] <- result_boxplot(x = result_matrices$rf_full,
               title = 'Prediction Performance',
               subtitle = 'Random Forest')

box_plots[['rf_small']] <- result_boxplot(x = result_matrices$rf_small,
               title = 'Prediction Performance',
               subtitle = 'Random Forest, small test sets')

box_plots[['svm_full']] <- result_boxplot(x = result_matrices$svm_full,
               title = 'Prediction Performance',
               subtitle = 'SVM')

box_plots[['svm_small']] <- result_boxplot(x = result_matrices$svm_small,
               title = 'Prediction Performance',
               subtitle = 'SVM, small test sets')

# Save figures
for (i in seq_along(box_plots)){
  file_name <- sprintf(
    'figures/assess_normalization/%s_grey.tiff',
    names(box_plots)[i]
  )
  ggsave(
    filename = file_name,
    plot = box_plots[[i]],
    width = 5,
    height = 5,
    dpi = 400
  )
}
# Save result metrics
for (i in seq_along(mean_mcc)){
  file_name <- sprintf(
    'results/summarized/%s.txt',
    names(mean_mcc)[i]
  )
  Encoding(x = mean_mcc[[1]]$output) <- 'UTF-8'
  fwrite(
    x = mean_mcc[[i]],
    file = file_name,
    sep = '\t',
    row.names = TRUE
  )
}
