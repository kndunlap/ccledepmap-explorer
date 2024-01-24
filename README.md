Hello! I'm Kyle and I'm here with ccledepmap-explorer. This is a way to visualize and interact with the dependency scores found in the Achilles DEPMAP dataset.
A lot of this work was inspired by an excellent paper out of the lab of Lydia Finley. They analyzed CRISPR essentiality scores and observed that genes involved in the TCA cycle gather into distinct clusters. As you can see by the figure below.
<img width="663" alt="Finley_Fig1A" src="https://github.com/kndunlap/ccledepmap-explorer/assets/61035909/024fc69b-20b8-474b-81e4-682ea59f0207">

I made it my mission to first try and re-create this figure. Once I had done that, I was planning on making a network graph for the Urea Cycle genes that I often work with.
The first thing that I had to do was collect some datasets. These are stored on the Utah Biochem Fileserver.
```
ccle_metadata <- read_csv("0 - sample_info.csv")
gene_effect <- read_csv("1 - CRISPRGeneEffect.csv")
gene_effect_metadata <- read_csv("2 - gene_effect_withinfo.csv")
gene_effect_cor_matrix <- read_csv("3 - gene_effect_cor_matrix.csv")
gene_dependency <- read_csv("4 - CRISPRGeneDependency.csv")
gene_dependency_metadata <- read_csv("5 - gene_dependency_withinfo.csv")
gene_dependency_cor_matrix <- read_csv("6 - gene_dependency_cor_matrix.csv")
```
