#!/bin/bash

# job script used to run cell segmentation on membranes (in this case, for P51, but completed on each patient's data)
source myconda
conda activate snakemakePipeline
python -m cellpose --dir /data/Zhaolab/1_AMLCosMx/Final_scripts/1_Normalization/0_NormalizedImg/P51_R1158_S1_Normalized/ --pretrained_model /gpfs/gsfs12/users/Zhaolab/1_AMLCosMx/v7_training/train_membrane_1_11_23/train/models/cellpose_residual_on_style_on_concatenation_off_train_2023_01_12_16_14_47.644146 --chan 2 --chan2 0 --diameter 0. --save_tif --use_gpu --savedir /data/Zhaolab/1_AMLCosMx/Final_scripts/2_Segmentation/1_MembraneSegmentation/P51_v7_output_membrane/prediction_model_1_11_23/ --verbose --no_npy
