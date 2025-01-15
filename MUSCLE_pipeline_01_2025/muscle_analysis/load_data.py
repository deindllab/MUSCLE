import tkinter.filedialog as fd
from skimage import transform
import numpy as np
import tkinter as tk
from tkinter import messagebox
import copy

from . import pks_from_img

from .constants import *


def poly_transf():
    """This function is used for uploading the polynomial transformation which defines the relation between Green and Red channels

    Args:
        base_path (string): The absolute path of the input files.
       
    Returns:
        tr_G2R: forward transformation (it describes transformation from green to red channel),
        tr_R2G: reverse transformation (it describes transformation from red to green channel)
    """
    tr_G2R = transform.PolynomialTransform()
    file_path = fd.askopenfilename(title = "Choose the forward transform file (Green to Red)")
    tr_G2R.params = np.load(file_path)
    tr_R2G = transform.PolynomialTransform()
    file_path = fd.askopenfilename(title = "Choose the inverse transform file (Red to Green)")
    tr_R2G.params = np.load(file_path) 
    return (tr_G2R, tr_R2G)


def extract_pos_info(POS, res = False):
    """
    This function is used for extraction the data from the position list. 
    If res = True, then it also allows to extract specific data from the position list.

    Args:
        POS: position list
        res: boolean
       
    Returns:
        labels - labels of the position
        posX - x coordinate of the position
        posY - y coordinate of the position 
                  OR
        labels_res - labels of the position, which has been chosen based on criteria
        posX_res -  x coordinate of the position, which has been chosen based on criteria
        posY_res -  y coordinate of the position, which has been chosen based on criteria
    """ 

    labels = [P['LABEL'] for P in POS]
    posX = [P['DEVICES'][1]['X'] for P in POS]
    posY = [P['DEVICES'][1]['Y'] for P in POS]
    pos = [posX,posY]
    if stage_rotation == 1:
        posx = stage_x_cor*np.array(pos[0])
        posy = stage_y_cor*np.array(pos[1])
    if stage_rotation == -1:
        posx = stage_x_cor*np.array(pos[1])
        posy = stage_y_cor*np.array(pos[0])
    if res:
        minx = min(posx)
        miny = min(posy)
        deltax, deltay = pks_from_img.get_translation()
        print(deltax,deltay)
        # print(upper_left_um_x, upper_left_um_y)
        labels_res = []
        posX_res = []
        posY_res = []
        for i,x in enumerate(labels):

            x = int(deltax-(x_border-converted_x)/2 + (posx[i] - minx)*1000/MiSeq_pixel)
            y = int(deltay-(y_border-converted_y)/2 + (posy[i] - miny)*1000/MiSeq_pixel)

            if (x > 0 ) and (y > 0) and (y+y_border< tile_size_y) and (x+x_border < tile_size_x):
                    labels_res.append (labels[i])
                    posX_res.append(x)
                    posY_res.append(y)
        #print (labels_res)
        print('Number of overlapping FOVs: ',len(labels_res))
        return (labels_res, posX_res, posY_res)

    else:
        return (labels, posx , posy)
