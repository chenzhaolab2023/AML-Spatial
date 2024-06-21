# Spatial Analysis: Prepare ligand-recptor input, summarize output

This directory contains the jupyter notebook "Prepare_LigandReceptor_input.ipynb" to quantify the number of each cell type at at distances of less than 5 microns or greater than 30 microns from each leukemia cell. The resulting file was then used as an input to run the pseudo-bulk ligand-receptor analysis (executed by Gege Gui in R). 

The ligand receptor analysis was perfomed in 2 directions: once with the leukemia cells assumed to express the ligands while the neighboring cells expressed the receptor and once where the leukemia were cells assumed to express the receptor while the neighboring cells expressed the ligand. The results were visualized using the notebook "Visualize_LRpairs.ipynb"
