# Cell Typing: ID and remove RBCs

This directory contains jupyter notebooks and scripts used to identify and remove RBCs. First, the snakefile is executed to crop single cell images and extract their features using EfficientNet. Then for each patient individually, at each timepoint individually, we cluster the features using Phenograph and identify the RBC cluster, overlaying the predicted RBCs on the original image.
