library(readr)
library(dplyr)
library(ggplot2)

busco_output <- list.files("output/busco",
                           recursive = TRUE,
                           pattern = ".tsv",
                           full.names = TRUE)
busco_fulltables <- grep("full_table", busco_output, value = TRUE)

busco_results_list <- lapply(busco_fulltables,
                             readr::read_tsv,
                             skip=4)

names(busco_results_list) <- gsub(".*full_table_(.+).tsv", "\\1", busco_fulltables)

busco_all <- dplyr::bind_rows(busco_results_list, .id = "filename")

plot_data <- busco_all %>%
  group_by(filename, Status) %>%
  summarise(n = n()) %>%
  group_by(filename) %>%
  mutate(percentage = n/sum(n)*100)

#order variables
status_order <- c("Complete", "Fragmented", "Duplicated", "Missing")
plot_data$Status <- factor(plot_data$Status, levels = rev(status_order))

ggplot(plot_data, aes(x = filename, y = percentage, fill = Status)) +
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_col()
