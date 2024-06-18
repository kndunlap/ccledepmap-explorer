# 1. Run these commands first - install the packages if you haven't ---------------------------------------------
# install.packages("tidyverse")
# install.packages("factoextra")
# install.packages("ggrepel")
library(tidyverse)
library(factoextra)
library(ggrepel)


# 2. Code - Run this but don't change anything ----------------------------

cluster_graph2 <- function(file, list, line_cutoff = 0.25, title) {
  
  gene_dependency_cor_matrix <- read_csv(file, col_types = cols(.default = "d", ...1 = "c"))|>
    rename(Gene = ...1)
  
  colnames(gene_dependency_cor_matrix) <- sub("\\s*\\(\\d+\\)", "", colnames(gene_dependency_cor_matrix))
  
  example <- read_csv(list)
  
  names(example)[1] <- 'Gene'
  
  list <- as.vector(unlist(example[,1]))
  
  list1 <- list[!(list %in% colnames(gene_dependency_cor_matrix))]
  print(list1)
  if(any(!list %in% colnames(gene_dependency_cor_matrix))) {
    stop(paste0("Sorry, I don't recognize the following gene(s):\n ", 
                list1, 
                " \nTry again! It might not exist in this dataset, or you have called it by the wrong name."))
  }
  
  cor_matrix <- gene_dependency_cor_matrix |>
    select(-Gene) |>
    cor() 
  
  cor_matrix_tibble <- as.tibble(cor_matrix)
  cor_matrix_final_tibble <- tibble(Gene = colnames(cor_matrix_tibble), cor_matrix_tibble)
  
  list1 <- list
  list1 <- sort(list1)
  
  tca_matrix <- cor_matrix_final_tibble |>
    select(all_of(list1), Gene) |>
    filter(Gene %in% list1) |>
    select(!Gene)
  
  data.pca <- princomp(tca_matrix)
  summary(data.pca)
  
  pca_data <- data.pca$loadings[, 1:2]
  pca_data
  
  pca_data <- as.data.frame(pca_data)
  pca_data
  
  gene_names <- rownames(pca_data)
  gene_names
  
  pca_data_complete <- cbind(pca_data[, 1:2], Gene = gene_names)
  pca_data_complete
  
  pca_data <- as.tibble(pca_data_complete)
  
  has_second_column <- ncol(example) > 1
  
  if (has_second_column) {
    pca_data <- pca_data |>
      left_join(example)
  }
  
  result_df <- tca_matrix |>
    as.data.frame() |>
    rownames_to_column(var = "Gene1") |>
    pivot_longer(
      cols = !1,
      names_to = "Gene2",
      values_to = "Value"
    ) |>
    select(!1)
  
  gene_df <- data.frame(Genes = rep(list1, each = length(list1)))
  gene_df <- as_tibble(gene_df)
  final_df <- as_tibble(cbind(gene_df, result_df))
  final_df <- final_df |> rename(Gene1 = Genes)
  
  final_df <- final_df |>
    mutate(plot_line_if = ifelse(Value > line_cutoff, Value, 0))
  
  if (has_second_column) {
    p <- pca_data |>
      ggplot(aes(x = Comp.1, y = Comp.2)) +
      geom_point(aes(size = 6, alpha = 0.8, color = Group)) +
      geom_text_repel(aes(label = Gene), size = 3.5, force = 11) +
      geom_segment(data = final_df, aes(x = pca_data[match(Gene1, pca_data$Gene), ]$Comp.1,
                                        y = pca_data[match(Gene1, pca_data$Gene), ]$Comp.2,
                                        xend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.1,
                                        yend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.2,
                                        alpha = plot_line_if)
      ) +
      scale_alpha_continuous(range = c(0, 1), guide = FALSE) +
      theme_void() +
      labs(title = paste0(title, " Networking Chart"),
           subtitle = paste0("Line cutoff r < ", line_cutoff, "\n Line Opaqueness Corresponds to Correlation Strength")) +
      theme(plot.title = element_text(hjust = 0.5, size = 20), 
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            legend.position = c(0.07, 0.98)) +
      guides(size = "none", 
             color = guide_legend(override.aes = list(size = 6)))
  } else {
    p <- pca_data |>
      ggplot(aes(x = Comp.1, y = Comp.2)) +
      geom_point(aes(size = 6, alpha = 0.8, color = "red")) +
      geom_text_repel(aes(label = Gene), size = 3.5, force = 11) +
      geom_segment(data = final_df, aes(x = pca_data[match(Gene1, pca_data$Gene), ]$Comp.1,
                                        y = pca_data[match(Gene1, pca_data$Gene), ]$Comp.2,
                                        xend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.1,
                                        yend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.2,
                                        alpha = plot_line_if)
      ) +
      scale_alpha_continuous(range = c(0, 1), guide = FALSE) +
      theme_void() +
      labs(title = paste0(title, " Networking Chart"),
           subtitle = paste0("Line cutoff r < ", line_cutoff, "\n Line Opaqueness Corresponds to Correlation Strength")) +
      theme(plot.title = element_text(hjust = 0.5, size = 20), 
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            legend.position = c(0.07, 0.98)) +
      guides(size = "none")
  }
  
  print(p)
}

# 4. Run command - give 4 inputs. ---------------------------------------


# a. file name (Must be a CSV in your directory)
# b. genelist and graph grouping if you want (Must be a CSV in your directory)
# c. A correlation cutoff. Default is 0.25
# d. A title for your graph, must be in quotes.

# Example function call is below 

cluster_graph2("CRISPRGeneDependency.csv", "example_less.csv", 0.30, "Urea Cycle ")

