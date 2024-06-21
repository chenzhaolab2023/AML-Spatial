#!/bin/bash

# job script used to run cell segmentation on nuclei (in this case, for P51, but completed on each patient's data)
source myconda
module load CUDA/11.3
module load cellpose
cellpose --dir /data/Zhaolab/1_AMLCosMx/Final_scripts/1_Normalization/0_NormalizedImg/P51_R1158_S1_Normalized/ --pretrained_model /data/Zhaolab/v7_training/nucleus_train_12_19/train/models/cellpose_residual_on_style_on_concatenation_off_train_2022_12_19_22_59_20.888505 --chan 3 --chan2 0 --diameter 0. --save_tif --use_gpu --savedir /data/Zhaolab/1_AMLCosMx/Final_scripts/2_Segmentation/0_NuclearSegmentation/P51_v7_output/labels_predicted/ --no_npy --verbose
