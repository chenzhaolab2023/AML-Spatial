#!/bin/bash

source myconda
conda activate snakemakePipeline
cd /data/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/4_CelltypesByTranscript/
snakemake -s snakefile_overlay_celltypes --profile /data/Zhaolab/1_AMLCosMx/v7_training/snakemake_profile/
