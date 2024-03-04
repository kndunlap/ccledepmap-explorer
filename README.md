# Overview
Hello! I'm Kyle and I'm here with ccledepmap-explorer. This README will instruct you on how to create a network plot using Dependency scores from CERES. This allows you to make plots similar to the one in Fig. 1A of this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8934290/. This tool uses Principal Component Analysis (PCA) to cluster genes together. Below is a re-creation of this figure in the link.

<img width="548" alt="Picture1" src="https://github.com/kndunlap/ccledepmap-explorer/assets/61035909/0661d973-d2cb-428a-b68a-e441cd1cd491">

I have two scripts available for now - one is if you already have the correlation matrix, and the other is if you only have the CERES score sheet.

To use the first script, you need to find "script1_cluster_graph.R" within this repository. For input in the command, you need to give R a correlation matrix sheet and a list of genes. If you want the one that I made from the most updated CERES Dependency Dataset, download the ~2.2 GB file from this link: https://drive.google.com/file/d/1EQsOu6X4tc_sNnqsphnR6tAxMy7SiijW/view?usp=drive_link. Then, that can be your input when you run the code.

To use the second script, "script2_cluster_graph.R", you can give R a dependency list of your choice. You also will give a list of genes. You might want to do this if you are wanting to customize the cell lines or experimental setup that you input. If you do that you'll probably have to manipulate the input sheet somehow but I can help with that. More than likely, you will find the dependency lists from this link: https://depmap.org/portal/download/all/


