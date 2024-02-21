# Overview
Hello! I'm Kyle and I'm here with ccledepmap-explorer. This README will instruct you on how to create a network plot using Dependency scores from CERES. This allows you to make plots similar to the one in Fig. 1A of this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8934290/

I have one script available for now - one is if you already have the correlation matrix, and the other is if you only have the CERES score sheet.

To use the first script, you need to find "script1_cluster_graph.R" within this repository. For input in the command, you need to give R a correlation matrix sheet and a list of genes. If you want the one that I made from the most updated CERES Dependency Dataset, download the ~2.2 GB file from this link: https://drive.google.com/file/d/1EQsOu6X4tc_sNnqsphnR6tAxMy7SiijW/view?usp=drive_link. Then, that can be your input when you run the code.

I'm still working on this part - but to use the second script, you can give R a dependency list of your choice. You might want to do this if you are only looking in a subset of cell lines for example. If you do that you'll probably have to manipulate the input sheet somehow but I can help with that.


