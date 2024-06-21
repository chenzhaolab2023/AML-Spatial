import numpy as np
import cv2
from skimage import measure
from shapely import Polygon, distance
import tifffile
import pandas as pd

## This script computes the shortest edge-to-edge distance between cells in each FOV

def prepare_box_for_contours(box, shape, pad=3):
    """Marginally pads a bounding box so that object boundaries
    are not on cropped image patch edges.
    """
    for i in range(2):
        box[i] = min(0, box[i] - pad)
        box[i+2] = max(shape[i], box[i] + pad)
        
    slices = tuple([slice(box[i], box[i+2]) for i in range(2)])
    top_left = np.array(box[:2])[None] # (1, 2)
    return slices, top_left

def make_polygons_from_mask(mask):
    """Constructs a polygon for each object in a mask. Returns
    a dict where each key is a label id and values are shapely polygons.
    """
    polygons = {}
    i = 0
    for rp in measure.regionprops(mask):
        # Faster to compute contours on small cell tiles than the whole image
        box_slices, box_top_left = prepare_box_for_contours(list(rp.bbox), mask.shape)
        label_mask = mask[box_slices] == rp.label
        #print(i)
        label_cnts = np.concatenate(
            measure.find_contours(label_mask), axis=0
        )

        polygons[rp.label] = Polygon(label_cnts + box_top_left)
        i += 1
    return polygons

def pairwise_polygon_distance(polygons_dict):
    """Computes pairwise distance between all polygons in
    a dictionary. Returns a dictionary of distances.
    """
    distances = {l: {} for l in polygons_dict.keys()}
    for i in polygons_dict.keys():
        for j in polygons_dict.keys():
            # nested loop is slow but we cache results
            # to eliminate duplicate work
            if i != j and distances[i].get(j) is None:
                distances[i][j] = distance(polygons_dict[i], polygons_dict[j])
                
    return distances

def pairwise_polygon_distance_2f(polygons_dict):
    """Computes pairwise distance between all polygons in
    a dictionary. Returns a dictionary of distances.
    """
    distances = {l: {} for l in polygons_dict.keys()}
    for i in polygons_dict.keys():
        for j in polygons_dict.keys():
            # nested loop is slow but we cache results
            # to eliminate duplicate work
            if i != j and distances[i].get(j) is None:
                distances[i][j] = format(distance(polygons_dict[i], polygons_dict[j]), '.2f')
                
    return distances

def get_nn_distance(key, distances_dict):
    """Returns the nearest neighbor for a polygon
    along with the distance.
    """
    return min(distances_dict[key].items(), key=lambda x: x[1])


# load mask file
mask = tifffile.imread(snakemake.input[0])

# zero-pad edges
mask_padded = np.pad(mask, (5, 5), 'constant')

# compute shortest distances between all cells
polygons_dict = make_polygons_from_mask(mask_padded)
distances_dict = pairwise_polygon_distance_2f(polygons_dict)

# save dataframe of shortest distances between all cells
min_dists = pd.DataFrame.from_dict(distances_dict, orient='index')
min_dists.to_csv(snakemake.output[0])