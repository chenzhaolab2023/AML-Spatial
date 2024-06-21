import numpy as np
import tifffile
import os
import csv
import cv2
import pandas

def generate_counts(nano_counts_tx, FOV, masks, img_raw_z1, img_raw_z2, img_raw_z3, img_raw_z4, img_raw_z5, img_raw_z6, img_raw_z7, img_raw_z8, img_raw_z9):
    
    # read in nanostring counts table
    counts = pandas.read_csv(nano_counts_tx, index_col=False)
    
    # input which fov (called fov7 here, but each individual FOV is read here)
    fov7 = counts[counts.fov == int(FOV)]
    
    # Revise columns
    fov7 = fov7.drop(columns = ['CellComp'])
    fov7['Size (pixels)'] = 0
    fov7['CD298_total'] = 0
    fov7['CD298_mean'] = 0
    fov7['CD298_median'] = 0
    fov7['CD3_total'] = 0
    fov7['CD3_mean'] = 0
    fov7['CD3_median'] = 0
    fov7['CD45_total'] = 0
    fov7['CD45_mean'] = 0
    fov7['CD45_median'] = 0
    fov7['CD34_total'] = 0
    fov7['CD34_mean'] = 0
    fov7['CD34_median'] = 0
    fov7['DAPI_total'] = 0
    fov7['DAPI_mean'] = 0
    fov7['DAPI_median'] = 0
    
    # remove original cell id label
    fov7.loc[:,'cell_ID'] = 0
    
    # load mask file
    labels = tifffile.imread(masks)
    
    # flip labels
    labels = np.flipud(labels)
    
    # assign new cell IDs
    for i in range(fov7.shape[0]):
        fov7.cell_ID.iloc[i] = labels[round(fov7.y_local_px.iloc[i]), round(fov7.x_local_px.iloc[i])]
    
    # Read in OME-TIF with raw values for all z slices
    img_raw_z1 = tifffile.imread(img_raw_z1)
    
    img_raw_z2 = tifffile.imread(img_raw_z2)
    
    img_raw_z3 = tifffile.imread(img_raw_z3)
    
    img_raw_z4 = tifffile.imread(img_raw_z4)
    
    img_raw_z5 = tifffile.imread(img_raw_z5)
    
    img_raw_z6 = tifffile.imread(img_raw_z6)
    
    img_raw_z7 = tifffile.imread(img_raw_z7)
    
    img_raw_z8 = tifffile.imread(img_raw_z8)
    
    img_raw_z9 = tifffile.imread(img_raw_z9)
    
    z_stack = np.stack([img_raw_z1, img_raw_z2, img_raw_z3, img_raw_z4, img_raw_z5, img_raw_z6, img_raw_z7, img_raw_z8, img_raw_z9], axis=0)
    
    # z-projection for each channel
    z_avg = np.sum(z_stack, axis=0, dtype=int) / z_stack.shape[0]
    
    # Record protein values from each channel
    CD298 = z_avg[0,:,:]
    CD298 = np.flipud(CD298)
    
    CD3 = z_avg[1,:,:]
    CD3 = np.flipud(CD3)

    CD45 = z_avg[2,:,:]
    CD45 = np.flipud(CD45)
    
    CD34 = z_avg[3,:,:]
    CD34 = np.flipud(CD34)
    
    DAPI = z_avg[4,:,:]
    DAPI = np.flipud(DAPI)
    
    
    # Add protein levels for each cell
    for i in range(1, fov7['cell_ID'].max()+1):
        fov7.loc[fov7['cell_ID']==i, 'CD298_total'] = sum(CD298[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD3_total'] = sum(CD3[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD45_total'] = sum(CD45[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD34_total'] = sum(CD34[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'DAPI_total'] = sum(DAPI[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD298_mean'] = np.average(CD298[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD3_mean'] = np.average(CD3[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD45_mean'] = np.average(CD45[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD34_mean'] = np.average(CD34[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'DAPI_mean'] = np.average(DAPI[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD298_median'] = np.median(CD298[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD3_median'] = np.median(CD3[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD45_median'] = np.median(CD45[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'CD34_median'] = np.median(CD34[labels == i])
        fov7.loc[fov7['cell_ID']==i, 'DAPI_median'] = np.median(DAPI[labels == i])
        
    # Add cell size
    for i in range(1, fov7['cell_ID'].max()+1):
        fov7.loc[fov7['cell_ID']==i, 'Size (pixels)'] = np.sum(labels == i)
        
    
    return fov7

    

    

    

    
