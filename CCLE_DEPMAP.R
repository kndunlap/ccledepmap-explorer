
# Load Required Packages --------------------------------------------------


library(tidyverse)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)


# Import Files ------------------------------------------------------------


ccle_metadata <- read_csv("0 - sample_info.csv")
gene_effect <- read_csv("1 - CRISPRGeneEffect.csv")
gene_effect_metadata <- read_csv("2 - gene_effect_withinfo.csv")
gene_effect_cor_matrix <- read_csv("3 - gene_effect_cor_matrix.csv")
gene_dependency <- read_csv("4 - CRISPRGeneDependency.csv")

gene_dependency_metadata <- read_csv("5 - gene_dependency_withinfo.csv")

gene_dependency_cor_matrix <- read_csv("6 - gene_dependency_cor_matrix.csv")
gene_dependency_cor_matrix <- gene_dependency_cor_matrix |> select(!1)


# Framework For Cleaning and Joining --------------------------------------

ccle_metadata <- as.tibble(ccle_metadata)
ccle_metadata
ccle_metadata <- ccle_metadata |>
  rename(ID = DepMap_ID)
ccle_metadata


gene_dependency <- as.tibble(gene_dependency)

names(gene_dependency) <- sub("\\.\\.\\d+\\.", "", names(gene_dependency))
gene_dependency <- gene_dependency |>
  rename(ID = X)

write.csv(gene_dependency, "4 - CRISPRGeneDependency.csv")

gene_dependency_metadata <- ccle_metadata |>
  left_join(gene_dependency) 

View(gene_dependency_metadata)

gene_dependency_metadata <- gene_dependency_metadata |>
  filter(!is.na(AARS1))

View(gene_dependency_metadata)

write.csv(gene_dependency_metadata, "5 - gene_dependency_withinfo.csv")


# Framework for Making Correlation Matrix ---------------------------------

gene_dependency_cormatrix <- gene_dependency_metadata |>
  select(!1:29) |>
  cor()

gene_dependency_cormatrix_tibble <- as.tibble(gene_dependency_cormatrix)

gene_dependency_cormatrix_goodnames <- tibble(Gene = colnames(gene_dependency_cormatrix_tibble), gene_dependency_cormatrix_tibble)

gene_dependency_cor_matrix <- gene_dependency_cormatrix_goodnames
write.csv(gene_dependency_final, "gene_dependency_cor_matrix.csv")

# Framework for Making Correlation Matrix within a certain cancer type -----

gene_dependency_cor_matrix_colon <- gene_dependency_metadata |>
  relocate(primary_disease) |>
  select(!2:29) |>
  filter(primary_disease == "Colon/Colorectal Cancer") |>
  select(!1) |>
  cor()

gene_dependency_cor_matrix_colon_tibble <- as.tibble(gene_dependency_cor_matrix_colon)

gene_dependency_cor_matrix_colon_goodnames <- tibble(Gene = colnames(gene_dependency_cor_matrix_colon_tibble), gene_dependency_cor_matrix_colon_tibble)

gene_dependency_cor_matrix_colon <- gene_dependency_cor_matrix_colon_goodnames
write.csv(gene_dependency_final, "gene_dependency_cor_matrix.csv")

# Functions ---------------------------------------------------------------


# Get a List of Correlation Between Two Genes in all Cancer Types (output: 32x2) ---------

gene_dependency_metadata |>
  select(SLC7A1, SLC7A5, primary_disease) |>
  group_by(primary_disease) |>
  summarize(
    correlation = cor(SLC7A1, SLC7A5)
  ) |>
  arrange(desc(correlation)) |>
  print(n = Inf)


# Find what gene is the least co-essential with your gene  -----------------

least_essential <- function(Gene_name, howmany = 20) {
  gene_dependency_cor_matrix |>
  select({{Gene_name}}, Gene) |>
  arrange({{Gene_name}}) |>
  print(n = howmany)
}

least_essential(SLC7A1, 55)


# # Find what gene is the most co-essential with your gene  ------ --------

most_essential <- function(Gene_name, howmany = 20) {
  gene_dependency_cor_matrix |>
    select({{Gene_name}}, Gene) |>
    arrange(desc({{Gene_name}})) |>
    print(n = howmany)
}

most_essential(MDH1)

# Get a matrix of correlations with genes of your choice. -------------

list1 <- c("ACACA", "ACLY", "ACO2", "AHCY", "ALDOA", "CS", "DHFR2", "DLAT", "DLD", "DLST", "ENO1", "FASN", "GAPDH", "GPI", "HK2", "IDH3A", "MAT2A", "MDH2", "MTHFD1", "OGDH", "PDHA1", "PDHB", "PFKP", "PGK1", "PKM", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SLC25A1")
tca_matrix <- gene_dependency_cor_matrix |>
  select(Gene, ACACA, ACLY, ACO2, AHCY, ALDOA, CS, DHFR2, DLAT, DLD, DLST, ENO1, FASN, GAPDH, GPI, HK2, IDH3A, MAT2A, MDH2, MTHFD1, OGDH, PDHA1, PDHB, PFKP, PGK1, PKM, SDHA, SDHAF2, SDHB, SDHC, SDHD, SLC25A1) |>
  filter(Gene %in% list1) |>
  select(!Gene)
  print(n = Inf)
 

# Check which cancers your gene is most/least essential in ----------------

most_essential_cancers <- function(Gene) {
  gene_dependency_metadata |>
  group_by(primary_disease) |>
  summarize(
    new = mean({{Gene}})
  ) |>
    arrange(desc(new)) |>
    print(n = Inf)
}
most_essential_cancers(SLC7A5)

# Check which genes have the highest mean essentiality and total st. dev --------

check_mean_stdev <- function(howmany) {
  gene_dependency_metadata |>
  select(!2:19) |>
  pivot_longer(
    cols = !(2),
    names_to = "gene",
    values_to = "essentiality"
  ) |>
  group_by(gene) |>
  summarize(
    mean = mean(essentiality),
    sd = sd(essentiality)
  ) |>
  arrange(desc(mean)) |>
    print(n = howmany)
}
check_mean_stdev(-50)
  

# Make a Correlation Matrix Graph -----------------------------------------------

ggcorrplot(tca_matrix)


# Re-Creation of Arnold, et al. Fig 1A ------------------------------------


list1 <- c("ACACA", "ACLY", "ACO2", "AHCY", "ALDOA", "CS", "DHFR2", "DLAT", "DLD", "DLST", "ENO1", "FASN", "GAPDH", "GPI", "HK2", "IDH3A", "MAT2A", "MDH2", "MTHFD1", "OGDH", "PDHA1", "PDHB", "PFKP", "PGK1", "PKM", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SLC25A1", "TYMS")
tca_matrix <- gene_dependency_cor_matrix |>
  select(Gene, ACACA, ACLY, ACO2, AHCY, ALDOA, CS, DHFR2, DLAT, DLD, DLST, ENO1, FASN, GAPDH, GPI, HK2, IDH3A, MAT2A, MDH2, MTHFD1, OGDH, PDHA1, PDHB, PFKP, PGK1, PKM, SDHA, SDHAF2, SDHB, SDHC, SDHD, SLC25A1, TYMS) |>
  filter(Gene %in% list1) |>
  select(!Gene)
print(n = Inf)

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
pca_data |>
  print(n = Inf)

descrpt_vect <- c("FA Metabolism", "FA Metabolism", "TCA Cycle", "1C Metabolism", "Glycolysis", "TCA Cycle", "1C Metabolism", "TCA Cycle",
                  "TCA Cycle", "TCA Cycle", "Glycolysis", "FA Metabolism", "Glycolysis", "Glycolysis", "Glycolysis", "TCA Cycle", "1C Metabolism",
                  "TCA Cycle", "1C Metabolism", "TCA Cycle", "TCA Cycle", "TCA Cycle", "Glycolysis", "Glycolysis", "Glycolysis", "TCA Cycle", "TCA Cycle",
                  "TCA Cycle", "TCA Cycle", "TCA Cycle", "FA Metabolism", "1C Metabolism")
pca_data$Pathway <- descrpt_vect


result_df <- tca_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene1") %>%
  gather(key = "Gene2", value = "Value", -Gene1) %>%
  mutate(Gene1 = as.numeric(Gene1),
         Gene2 = as.numeric(gsub("X", "", Gene2)))

genes_list <- c("ACACA", "ACLY", "ACO2", "AHCY", "ALDOA", "CS", "DHFR2", "DLAT", "DLD", "DLST", "ENO1", "FASN", "GAPDH", "GPI", "HK2", "IDH3A", "MAT2A", "MDH2", "MTHFD1", "OGDH", "PDHA1", "PDHB", "PFKP", "PGK1", "PKM", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SLC25A1", "TYMS")
gene_df <- data.frame(Gene1 = rep(genes_list, each = 32))
gene_df$Gene2 <- genes_list

gene_df <- as.tibble(gene_df)
result_df <- as.tibble(result_df)

final_df <- cbind(gene_df, result_df)
final_df <- final_df[,-3:-4]
final_df <- as.tibble(final_df)


final_df <- final_df |>
  mutate(Valuee = ifelse(Value > 0.25, Value, 0))

pca_data |>
  ggplot(aes(x = Comp.1, y = Comp.2)) +
  geom_point(size = 6, aes(color = Pathway), alpha = 0.6) +
  geom_text_repel(aes(label = Gene), fontface = "bold", size = 3) +
  geom_segment(data = final_df, aes(x = pca_data[match(Gene1, pca_data$Gene), ]$Comp.1,
                                    y = pca_data[match(Gene1, pca_data$Gene), ]$Comp.2,
                                    xend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.1,
                                    yend = pca_data[match(Gene2, pca_data$Gene), ]$Comp.2,
                                    color = "Lines", alpha = Valuee)
               ) +
  scale_color_manual(values = c("Glycolysis" = "#eccd00", "FA Metabolism" = "blue", "TCA Cycle" = "red", "1C Metabolism" = "darkgreen")) +
  theme(legend.position = c(0.75, 0.7)) +
  scale_alpha_continuous(range = c(0, 1), guide = FALSE) +
  theme_void() +
  labs(title = "TCA Cycle Clustering Chart",
                    subtitle = ("Line cutoff r < .25\nOpaqueness Corresponds to Correlation Strength")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 12))



# cluster_graph Function, output is a network diagram --------------------------------------------------------

uc_list <- c("ASL",	"ASS1",	"ARG1",	"OTC", "FH",	"CPS1",	"CEBPA",	"ARG2",	"NAGS",	"AGMAT",	"SLC25A2",	"SLC25A15",	"SLC7A1",	"SLC7A5",	"CAD",	"ASNS",	"NOS1",	"NOS2",	"NOS3")

slc7list <- c("SLC7A1", "SLC7A2", "SLC7A3", "SLC7A4", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC7A9", "SLC7A10", "SLC7A11", "SLC7A13")

tca_list <- c("ACACA", "ACLY", "ACO2", "AHCY", "ALDOA", "CS", "DHFR2", "DLAT", "DLD", "DLST", "ENO1", "FASN", "GAPDH", "GPI", "HK2", "IDH3A", "MAT2A", "MDH2", "MTHFD1", "OGDH", "PDHA1", "PDHB", "PFKP", "PGK1", "PKM", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SLC25A1", "TYMS")

cluster_graph <- function(list, line_cutoff = 0.25, title = "Random") {
list1 <- list
list1 <- sort(list1)
tca_matrix <- gene_dependency_cor_matrix |>
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
  geom_text_repel(aes(label = Gene), fontface = "bold", size = 4, force = 23) +
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

cluster_graph(tca_list, 0.35, "Urea Cycle")
