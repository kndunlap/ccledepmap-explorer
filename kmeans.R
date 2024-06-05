library(tidyverse)
library(factoextra)
library(ggrepel)

gene_dependency_cor_matrix <- read_csv("7 - gene_dependency_cor_matrix_rounded.csv", col_types = cols(.default = "d", Gene = "c"))

gene_dependency_cor_matrix <- gene_dependency_cor_matrix |> select(!1)

list <- c("ACACA", "ACLY", "ACO2", "AHCY", "ALDOA", "CS", "DHFR2", "DLAT", "DLD", "DLST", "ENO1", "FASN", "GAPDH", "GPI", "HK2", "IDH3A", "MAT2A", 
           "MDH2", "MTHFD1", "OGDH", "PDHA1", "PDHB", "PFKP", "PGK1", "PKM", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SLC25A1")

list1 <- list
list1 <- sort(list1)

tca_matrix <- gene_dependency_cor_matrix |>
  select(all_of(list1), Gene) |>
  filter(Gene %in% list1) |>
  select(!Gene)


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


wide_data <- final_df %>%
  pivot_wider(names_from = Gene2, values_from = Value)

wide_data_df <- as.data.frame(wide_data)

rownames(wide_data_df) <- wide_data_df$Gene1

wide_data_df <- wide_data_df %>% select(-Gene1)

result_matrix <- as.matrix(wide_data_df)

print(result_matrix)

km.out <- kmeans(result_matrix, centers = 2, nstart = 20)

data_with_clusters <- as.data.frame(result_matrix)
data_with_clusters$Cluster <- as.factor(km.out$cluster)
data_with_clusters

data1 <- km.out$centers
View(data1)

km.out$cluster
