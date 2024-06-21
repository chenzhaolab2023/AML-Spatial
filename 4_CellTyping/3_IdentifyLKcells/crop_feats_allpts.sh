#!/bin/bash

module load python/3.9
cd /data/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/
snakemake -s snakefile_crop_feats_CD34_allpts --profile /data/Zhaolab/1_AMLCosMx/v7_training/snakemake_profile/
