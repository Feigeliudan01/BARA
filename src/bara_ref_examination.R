# Robin Gradin
#
# Examine the performance of BARA as the number
# of reference samples are varied.

library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(reshape2)
library(stringr)

# Read the result files
result_files <- list.files(path = 'results/bara_ref_samples/',
                           pattern = '.RDS$',
                           full.names = TRUE)
results <- list()
writeLines(sprintf('%s result files', length(result_files)))
for (i in seq_along(result_files)){
  file_name <- str_split(string = basename(result_files[i]),
                         pattern = '_')
  model <- paste(file_name[[1]][1], '.', file_name[[1]][2], sep = '')
  n_ref <- file_name[[1]][3] %>% gsub(pattern = '.RDS$', replacement = '')
  results[[model]][[n_ref]] <- readRDS(
    file = result_files[i]
  )
}

# Summarize the performances
results_mean <- map(
  results,
  function(model){
    map(
      model,
      function(n_ref){
        map_dbl(
          n_ref,
          function(training_set){
            median(map_dbl(training_set, ~mean(unlist(.))))
          }
        )
      }
    )
  }
)
performances <- map(results_mean, function(model){
  df <- as.data.frame(do.call(rbind, model))
  df$n_ref <- rownames(df)
  df
})

# Prepare for visualization
perf_full <- rbind.data.frame(
  performances$knn.full,
  performances$rf.full,
  performances$svm.full
)
perf_full$model <- c(
  rep('kNN', 6),
  rep('Random Forest', 6),
  rep('SVM', 6)
)
perf_full <- melt(perf_full, id.vars = c('n_ref', 'model')) %>%
  group_by(n_ref, model) %>%
  summarize(
    'mean_mcc' = mean(value),
    'sd_mcc' = sd(value)
  )

# Create plot
perf_plot <- ggplot(perf_full, aes(x = n_ref, y = mean_mcc)) +
  geom_line(aes(colour = model, group = model)) +
  geom_point(aes(colour = model)) +
  geom_errorbar(
    mapping = aes(
      colour = model,
      ymin = mean_mcc - sd_mcc,
      ymax = mean_mcc + sd_mcc
    ),
    width = 0.5
  ) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_colour_grey() +
  theme(
    text = element_text(family = 'serif'),
    panel.background = element_blank(),
    axis.line = element_line(),
    legend.key = element_blank()
  ) +
  labs(
    x = 'Number of Reference Samples',
    y = 'Mathews Correlation Coefficient',
    colour = 'Prediction Model'
  )
  # Create nicer table
  summarized_table <- perf_full[order(perf_full$model, perf_full$n_ref), ]
  summarized_table$mean_mcc <- round(summarized_table$mean_mcc, digits = 2)
  summarized_table$sd_mcc <- round(summarized_table$sd_mcc, digits = 2)

# Save plot and data
ggsave(
  filename = 'figures/bara/ref_full.tiff',
  width = 6,
  height = 4,
  dpi = 400
)
fwrite(
  x = perf_full,
  file = 'results/summarized/bara_n_ref_samples.txt',
  sep = '\t'
)
fwrite(
  x = summarized_table,
  file = 'results/summarized/bara_n_ref_samples_rounded.txt',
  sep = '\t'
)