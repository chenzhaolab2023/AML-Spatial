# Cell Typing: ID leukemia cells

This directory contains jupyter notebooks and scripts used to identify leukemia cells. First, the snakefile is executed to crop single cell images and extract features from the CD34 channel using EfficientNet. Then for each patient individually, at each timepoint individually, we cluster the features using Phenograph and identify the leukemia cluster(s), overlaying the predicted leukemia cells on the original image.

P51 was processed slightly differently, with an added step to remove the background CD34 staining since there was a higher level in background CD34 stain than was seen in the other patients.
