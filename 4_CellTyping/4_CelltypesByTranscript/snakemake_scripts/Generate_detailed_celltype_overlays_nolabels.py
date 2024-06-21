import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import cv2
from PIL import Image
from PIL import ImageDraw
from skimage import data, color, io, img_as_float
from os import listdir
from os.path import isfile, join
import skimage.io
import skimage.util
from skimage import data, segmentation, util, measure
from skimage.segmentation import expand_labels

# create dictionary of patient names
name_dirs = {'p1':'P51', 'p2':'P52', 'p3':'P53', 'p4':'P56', 'p5':'P57', 'p6':'P58'}
toGegeID = {'P51':'p1', 'P52':'p2', 'P53':'p3', 'P56':'p4', 'P57':'p5', 'P58':'p6'}

# Load prediction files
pred_location = 'FromGege/celltype_v8-EP.csv'

# read in celltype predictions
df = pd.read_csv(pred_location)

# reformat columns
df3 = pd.concat([df['x'], df['Unnamed: 0'].str.split('_', expand=True)], axis=1)
df3 = df3.rename(columns={"x": "celltype_detail", 0: "PtID", 3: "cell_id"})
df3['FOV'] = df3[1].apply(lambda x: x[-2:])
df3['patients'] = df3['PtID'].apply(lambda x: toGegeID[x])
df3 = df3.drop(columns=[1,2])
morph_predicted = df3[['patients', 'celltype_detail', 'FOV', 'cell_id']]

# define types to remove before further analysis
types_to_remove = ['NoCellAssigned','RBC','SmallCell','Unknown']

# filter out NA cells
morph_predicted = morph_predicted[~morph_predicted['celltype_detail'].isin(types_to_remove)]

# get relevant cell type names
class_column = 'celltype_detail'
cluster_labels = list(set(morph_predicted[class_column].tolist()))

# read csv with RGB values for each cell type
colors = pd.read_csv('color_v8.csv', index_col=0)
colors = colors.T

# extract colors, put in BGR order
cluster_colors = {}
for celltype in cluster_labels:
    cluster_colors[celltype] = [colors['blue'].loc[celltype], colors['green'].loc[celltype], colors['red'].loc[celltype]]

# load patient and FOV number from snakefile
patient = snakemake.params[1]
fov = snakemake.params[0]


## create overlay

# subset for patient
one_pt = morph_predicted[morph_predicted['patients'] == patient]

# subset FOV        
one_fov = one_pt[one_pt['FOV'] == fov]

# Load masks
labels = tifffile.imread(snakemake.input[0])

blank = np.zeros([labels.shape[0], labels.shape[1], 3])
blank = blank.astype('uint8')

# get all mask values
fov_cells = np.unique(labels)

combo = one_fov.set_index('cell_id')

for value in cluster_labels:

    # create list of this cluster's cells
    pos_cells = combo[combo['celltype_detail'] == value].index.tolist()
    pos_cells = [int(i) for i in pos_cells]

    # set all labels of other clusters to zero
    labels_edited = labels.copy()
    for cell in fov_cells:
        if cell in pos_cells:
            pass
        else:
            labels_edited[labels_edited == cell] = 0

    # color cells by type
    blank[labels_edited>0] = cluster_colors[value]

    # add colored outline
    boundaries = segmentation.find_boundaries(labels_edited, connectivity=1, mode='inner', background=0)
    boundaries = boundaries.astype(int)
    blank[boundaries>0] = [0,0,0] # add black borders to cells


# save image
cv2.imwrite(snakemake.output[0], blank)
