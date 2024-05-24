ccle <- read_csv("5 - gene_dependency_withinfo.csv")

ccle1 <- ccle |>
  select(!1)

View(ccle1)

ccle1 |>
  group_by(primary_disease) |>
  summarize(
    corvalue = cor(SLC7A5, ASS1), 
    n = n()
  ) |>
  arrange(desc(corvalue)) |>
  print(n = Inf)

ccleprad <- ccle1 |>
  filter(primary_disease == "Prostate Cancer")
  
cor_value <- cor(ccleprad$ASS1, ccleprad$SLC7A5)
  
ggplot(ccleprad, aes(x = SLC7A5, y = ASS1)) + 
  geom_point(alpha = 0.2) +
  annotate("text", x = min(ccleprad$ASS1), y = max(ccleprad$SLC7A5),
            label = paste("PRAD", "Correlation Coefficient: r =", round(cor_value, 3)),
            hjust = 0, vjust = 0, size = 5.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) 
}


ccle1 |>
  group_by(primary_disease) |>
  count() |>
  arrange(n) |>
  print(n = Inf)
