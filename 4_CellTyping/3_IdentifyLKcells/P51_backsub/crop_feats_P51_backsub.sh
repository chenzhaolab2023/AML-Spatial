#!/bin/bash

module load python/3.9
cd /data/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/P51_backsub/
snakemake -s snakefile_crop_feats_CD34_P51_backsub --profile /data/Zhaolab/1_AMLCosMx/v7_training/snakemake_profile/
