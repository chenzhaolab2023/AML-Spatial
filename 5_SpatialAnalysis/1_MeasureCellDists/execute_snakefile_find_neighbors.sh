#!/bin/bash

# This runs the snakefile which submits cell edge-edge distance calculations to the cluster for each FOV independently to parallelize and speed up the process.
source myconda
conda activate snakemakePipeline
cd /data/Zhaolab/1_AMLCosMx/Final_scripts/5_SpatialAnalysis/1_MeasureCellDists/
snakemake -s snakefile_find_neighbors --profile /data/Zhaolab/1_AMLCosMx/v7_training/snakemake_profile/
