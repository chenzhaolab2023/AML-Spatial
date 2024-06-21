#!/bin/bash

module load python/3.9
cd /data/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/1_IdentifyRBCs_2channel/
snakemake -s snakefile_crop_feats_RGB_allpts --profile /data/Zhaolab/1_AMLCosMx/v7_training/snakemake_profile/
