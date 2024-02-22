# 1. Run these commands first ---------------------------------------------
library(tidyverse)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)


# 2. Define list of genes (each gene must be in quotes) ----------------------------------------------

uc_list <- c("ASL",	"ASS1",	"ARG1",	"OTC", "FH",	"CPS1",	"CEBPA",	"ARG2",	"NAGS",	"AGMAT",	"SLC25A2",	
             "SLC25A15",	"SLC7A1",	"SLC7A5",	"CAD", "ASNS",	"NOS1",	"NOS2",	"NOS3")


# 3. Code - Run this but don't change anything ----------------------------

cluster_graph <- function(file, list, line_cutoff = 0.25, title) {
  
  gene_dependency_cor_matrix <- read_csv(file, col_types = cols(.default = "d", ...1 = "c"))|>
    rename(Gene = ...1)
  
  cor_matrix <- gene_dependency_cor_matrix |>
    select(-Gene) |>
    cor() 
    
  colnames(cor_matrix) <- sub("\\s*\\(\\d+\\)", "", colnames(cor_matrix))
  
  cor_matrix_tibble <- as.tibble(cor_matrix)
  cor_matrix_final_tibble <- tibble(Gene = colnames(cor_matrix_tibble), cor_matrix_tibble)
  
  if(any(!list %in% colnames(gene_dependency_cor_matrix))) {
    stop("Sorry, I don't recognize one of your genes. Try again! It might not exist in this dataset, or you have called it by the wrong name.")
  }
  
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
  gene_df <- as.tibble(gene_df)
  final_df <- as.tibble(cbind(gene_df, result_df))
  final_df <- final_df |> rename(Gene1 = Genes)
  
  
  final_df <- final_df |>
    mutate(plot_line_if = ifelse(Value > line_cutoff, Value, 0))
  
  pca_data |>
    ggplot(aes(x = Comp.1, y = Comp.2)) +
    geom_point(size = 6, alpha = 0.4, color = "red") +
    geom_text_repel(aes(label = Gene), fontface = "bold", size = 4, force = 11) +
    geom_segment(data = final_df, aes(x = pca_data[match(Gene1, pca_data$Gene), ]$Comp.1,
                                      y = pca_data[match(Gene1, pca_data$Gene), ]$Comp.2,
                                      xend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.1,
                                      yend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.2,
                                      alpha = plot_line_if)
    ) +
    theme(legend.position = c(0.85, 0.79)) +
    scale_alpha_continuous(range = c(0, 1), guide = FALSE) +
    theme_void() +
    labs(title = paste0(title, " Networking Chart"),
         subtitle = paste0("Line cutoff r < ", line_cutoff, "\n Line Opaqueness Corresponds to Correlation Strength")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 12), legend.position = "none") 
  
}

# 4. Run command - give 4 inputs. ---------------------------------------


# a. file name(Must be a CSV in your directory)
# b. the name of the list, defined in Step 2
# c. A line cutoff. Default is 0.25
# d. A title for your graph, must be in quotes.

# Example function call is below 

cluster_graph("CRISPRGeneDependency.csv", uc_list, 0.25, "Urea Cycle")

