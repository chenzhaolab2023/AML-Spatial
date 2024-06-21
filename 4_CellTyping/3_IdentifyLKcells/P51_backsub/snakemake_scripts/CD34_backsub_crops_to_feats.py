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

def Generate_v2_singlechannel(Original_loc, Mask_loc, patient, fov, max_dim): # This function crops single cells segmented in mask image, putting them on a black background (each as own max_dim x max_dim x 3 numpy array)
    
    # NOTE: since in this case the input is a grayscale PNG of the background-subtracted CD34 channel, we triplicate it to create a psuedo-RGB
    
    ## inputs:
    # Original_loc = path to background-subtracted CD34 PNG image to be cropped
    # Mask_loc = path to cell segmentation mask for corresponding image
    # patient = patient ID
    # fov = FOV #
    # max_dim = width/height of the cropped image - should be greater than diameter of largest cell
    
    ## ouputs:
    # cell_names = list of names of every cell in cell_stack in order
    # cell_stack = stack of numpy arrays, each a crop of a single cell

    mask = cv2.imread(Mask_loc, 2)
    image = cv2.imread(Original_loc, 2)
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
            new_image = np.zeros([max_dim, max_dim])
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
                
            # states that if cell is too large for max_dimxmax_dim box, add it to large cell list and skip the crop    
            if cell_image[y0:y1, x0:x1].shape[0] > new_image[yo:yi, xo:xi].shape[0] or cell_image[y0:y1, x0:x1].shape[1] > new_image[yo:yi, xo:xi].shape[1]:
                large_cells.append(name+1)
            else:
                new_image[yo:yi, xo:xi] = cell_image[y0:y1, x0:x1]
                resized_img = cv2.resize(new_image, (224, 224)) # resize for efficientNet
                img_bgr = np.stack([resized_img, resized_img, resized_img], axis=2) # triplicate for efficientNet
                all_cells.append(img_bgr) # create list of all arrays, then stack them
                cell_names.append(patient + '_FOV' + fov + '_cell_' + str(name+1)) # create list of all names

        name += 1
        
    print(str(fov) + " large cells: " + str(large_cells))
    
    # stack for efficientNet input
    cell_stack = np.stack(all_cells, axis=0)

    return cell_names, cell_stack


def extract_features(img_stack): # saves 672 image features from B0 layer of efficient net for each single cell crop
    
    model = EfficientNetB0(weights='imagenet')
    
    x = preprocess_input(img_stack)

    preds=model.predict(x)

    layer_name = 'block6a_activation'

    intermediate_layer_model = Model(inputs=model.input, outputs=model.get_layer(layer_name).output)
    intermediate_output = intermediate_layer_model.predict(x)

    feats = intermediate_output

    while len(feats.shape) > 2:  # 2D mean spatial pooling
        feats = np.mean(feats, axis=1)

    return feats ## returns 672 element feature array for each cell

## code to execute

# generate single cell crops for all cells
cell_names, cell_stack = Generate_v2_singlechannel(snakemake.input[0], snakemake.input[1], snakemake.params[0], snakemake.params[1], 300)

# run feature extraction on each channel, save as DF
ch0_features = extract_features(cell_stack)
DF0 = pd.DataFrame(ch0_features, index=cell_names)

# save features to csv
DF0.to_csv(snakemake.output[0])


