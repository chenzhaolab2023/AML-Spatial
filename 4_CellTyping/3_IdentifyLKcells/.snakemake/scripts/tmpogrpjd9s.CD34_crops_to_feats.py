
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/usr/local/Anaconda/envs/py3.9/lib/python3.9/site-packages', '/home/annemb/.cache/snakemake/snakemake/source-cache/runtime-cache/tmp9x_akgo8/file/gpfs/gsfs12/users/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/snakemake_scripts', '/gpfs/gsfs12/users/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/snakemake_scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95v\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\xa0/data/Zhaolab/1_AMLCosMx/Final_scripts/1_Normalization/0_NormalizedImg/P57_R1158_S3_Normalized_DAPI_B2M_CD34/20220228_110038_S3_C902_P99_N99_F015_normalized.png\x94\x8c\xa6/data/Zhaolab/1_AMLCosMx/Final_scripts/2_Segmentation/3_NucMemMerging/P57_hybrid/labels_predicted_2_15_24/20220228_110038_S3_C902_P99_N99_F015_normalized_cp_masks.tif\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x11\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x17)}\x94\x8c\x05_name\x94h\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cx/data/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/P57/Normalized_CD34_feats/P57_FOV15_CD34_feats.csv\x94a}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c\x03P57\x94\x8c\x0215\x94e}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94h6a}\x94(h\r}\x94\x8c\x04fovs\x94K\x00N\x86\x94sh\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bhHh6ub\x8c\x07threads\x94K\x02\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x02K\x01J\xa0\x86\x01\x00J\x88t\x01\x00M\xe8\x03M\xba\x03\x8c\x16/lscratch/20453948/tmp\x94K<K\x01\x8c\x05v100x\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06mem_mb\x94K\x02N\x86\x94\x8c\x07mem_mib\x94K\x03N\x86\x94\x8c\x07disk_mb\x94K\x04N\x86\x94\x8c\x08disk_mib\x94K\x05N\x86\x94\x8c\x06tmpdir\x94K\x06N\x86\x94\x8c\x07runtime\x94K\x07N\x86\x94\x8c\x03gpu\x94K\x08N\x86\x94\x8c\tgpu_model\x94K\tN\x86\x94uh\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bh\\K\x02h^K\x01h`J\xa0\x86\x01\x00hbJ\x88t\x01\x00hdM\xe8\x03hfM\xba\x03hhhXhjK<hlK\x01hnhYub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x11sNt\x94bh\x12h\x15h\x17\x85\x94R\x94(h\x17)}\x94h\x1bh\x12sNt\x94bub\x8c\x06config\x94}\x94(\x8c\rrestart-times\x94K\x00\x8c\tjobscript\x94\x8c\x12slurm_jobscript.sh\x94\x8c\x07cluster\x94\x8c\x0cbw_submit.py\x94\x8c\x0ecluster-status\x94\x8c\x0cbw_status.py\x94\x8c\x13max-jobs-per-second\x94K\x01\x8c\x1cmax-status-checks-per-second\x94G?\x84z\xe1G\xae\x14{\x8c\x0blocal-cores\x94K\x02\x8c\x0clatency-wait\x94K\xf0\x8c\x04jobs\x94K\x19\x8c\nkeep-going\x94\x88\x8c\x10rerun-incomplete\x94\x88\x8c\x07verbose\x94\x88u\x8c\x04rule\x94\x8c\tfeats_P57\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cd/gpfs/gsfs12/users/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/snakemake_scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/gpfs/gsfs12/users/Zhaolab/1_AMLCosMx/Final_scripts/4_CellTyping/3_IdentifyLKcells/snakemake_scripts/CD34_crops_to_feats.py';
######## snakemake preamble end #########
import math
import numpy as np
import cv2
import os
import tifffile
import pandas as pd
import csv
from os import listdir
from os.path import isfile, join
from tensorflow.keras.applications import EfficientNetB0
from tensorflow.keras.preprocessing import image
from tensorflow.keras.applications.imagenet_utils import decode_predictions
from tensorflow.keras.applications.imagenet_utils import preprocess_input
from keras.models import Model

## define functions

def Generate_v2(Original_loc, Mask_loc, patient, fov, max_dim): 

    mask = cv2.imread(Mask_loc, 2)
    image = cv2.imread(Original_loc)
    x_length = np.shape(mask)[1]
    y_length = np.shape(mask)[0]
    num = int(np.max(mask))
    cell = []
    for i in range(0, num):
        cell.append([])
    for x in range(0, x_length - 1):
        for y in range(0, y_length - 1):
            if int(mask[y][x]) != 0:
                cell[int(mask[y][x]) - 1].append([x, y])
                
    all_cells = []
    cell_names = []
    
    large_cells = []
    name = 0
    for n in cell:
        if len(n) != 0:
            y0, x0, y1, x1 = np.amin(n, axis=0)[1], np.amin(n, axis=0)[0], np.amax(n, axis=0)[1], np.amax(n, axis=0)[0]
            cell_t = np.zeros_like(image)
            for xy in n:
                cell_t[xy[1]][xy[0]] = 1
            cell_image = image * cell_t
            new_image = np.zeros([max_dim, max_dim, 3])
            xo = int(max_dim/2 - (x1 - x0) * 0.5)
            yo = int(max_dim/2 - (y1 - y0) * 0.5)
            xi = max_dim - xo
            yi = max_dim - yo
            if (x1 - x0) % 2 != 0:
                xo = math.floor(max_dim/2 - (x1 - x0) * 0.5) + 1
                xi = max_dim+1 - xo
            if (y1 - y0) % 2 != 0:
                yo = math.floor(max_dim/2 - (y1 - y0) * 0.5) + 1
                yi = max_dim+1 - yo
                
            # states that if cell is too large for max_dim x max_dim box, add it to large cell list and skip the crop    
            if cell_image[y0:y1, x0:x1, :].shape[0] > new_image[yo:yi, xo:xi, :].shape[0] or cell_image[y0:y1, x0:x1, :].shape[1] > new_image[yo:yi, xo:xi, :].shape[1]:
                large_cells.append(name+1)
            else:
                new_image[yo:yi, xo:xi, :] = cell_image[y0:y1, x0:x1, :]
                
                # resize for efficientNet
                resized_img_ch0 = cv2.resize(new_image[:,:,0], (224, 224))
                resized_img_ch1 = cv2.resize(new_image[:,:,1], (224, 224))
                resized_img_ch2 = cv2.resize(new_image[:,:,2], (224, 224))
                
                # stack back to rgb (pseudo-rgb, red channel stacked x3)
                img_rrr = np.stack([resized_img_ch2, resized_img_ch2, resized_img_ch2], axis=2)
                all_cells.append(img_rrr) # create list of all arrays, then stack them
                cell_names.append(patient + '_FOV' + fov + '_cell_' + str(name+1)) # create list of all names

        name += 1
        
    print(str(fov) + " large cells: " + str(large_cells))
    
    # stack for efficientNet input
    cell_stack = np.stack(all_cells, axis=0)

    return cell_names, cell_stack


def extract_features(img_stack):
    
    model = EfficientNetB0(weights='imagenet')
    
    x = preprocess_input(img_stack)

    preds=model.predict(x)

    layer_name = 'block6a_activation'

    intermediate_layer_model = Model(inputs=model.input, outputs=model.get_layer(layer_name).output)
    intermediate_output = intermediate_layer_model.predict(x)

    feats = intermediate_output

    while len(feats.shape) > 2:  # 2D mean spatial pooling
        feats = np.mean(feats, axis=1)

    return feats

## code to execute

# generate single cell crops for all cells
cell_names, cell_stack = Generate_v2(snakemake.input[0], snakemake.input[1], snakemake.params[0], snakemake.params[1], 300)

# run feature extraction on each channel, save as DF
ch0_features = extract_features(cell_stack)
DF0 = pd.DataFrame(ch0_features, index=cell_names)

# save features to csv
DF0.to_csv(snakemake.output[0])


