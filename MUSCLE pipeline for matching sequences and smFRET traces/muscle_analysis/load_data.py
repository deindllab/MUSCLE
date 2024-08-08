import tkinter.filedialog as fd
from skimage import transform
import numpy as np
import tkinter as tk
from tkinter import messagebox

from . import pks_from_img


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

def apriori_transf(x_border,y_border):
    """This function is used for creating the apriori transformation and it's further usage in code

    Args:
        x_border,y_border (int): 
       
    Returns:
        apriori_tr_original - original apriori transformation
        apriori_tr - centred transformation 
        apriori_tr_inv - inversed and centred transformation 
    """    
    
    FRET_pixel = 204 # nm/pixel
    Illumina_pixel = 340 #nm/pixel 
    scale = FRET_pixel/ Illumina_pixel
    apriori_tr_original = transform.SimilarityTransform(scale = scale)
    
    apriori_tr = transform.SimilarityTransform(scale = scale)
   # apriori_tr_original = transform.SimilarityTransform()
   # apriori_tr_original.params = np.array([[ 5.98278480e-01, -2.04811207e-03,  0],
   #  [ 2.04811207e-03,  5.98278480e-01,  0],
   #  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
   # apriori_tr = transform.SimilarityTransform()
   # apriori_tr.params = np.array([[ 5.98278480e-01, -2.04811207e-03,  0],
   #  [ 2.04811207e-03,  5.98278480e-01,  0],
   #  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    # Centering the tranformation to the x_border and y_border
    delta = 0.5*np.subtract([x_border,y_border],apriori_tr([512,256]))
    dx, dy = delta[0]
    apriori_tr.params[0,2] = dx
    apriori_tr.params[1,2] = dy
    apriori_tr_inv = transform.SimilarityTransform(scale = 1/scale)
   # apriori_tr_inv = transform.SimilarityTransform()
   # apriori_tr_inv.params = np.array([[ 1.67127284e+00,  5.72133913e-03,  0],
# [-5.72133913e-03,  1.67127284e+00,  0],
 #[ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

     # Centering the tranformation to the x_border and y_border
    delta = 0.5*np.subtract([512,256],apriori_tr_inv([x_border,y_border]))
    dx, dy = delta[0]
    apriori_tr_inv.params[0,2] = dx
    apriori_tr_inv.params[1,2] = dy
    
    
    return (apriori_tr_original, apriori_tr, apriori_tr_inv)

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
    
    x_border = 500
    y_border = 300
    labels = [P['LABEL'] for P in POS]
    posX = [P['DEVICES'][1]['X'] for P in POS]
    posY = [P['DEVICES'][1]['Y'] for P in POS]
    if res:
        #coords_t_seq = np.array([-2988.2,205.6])
        
        
        
        upper_left_um_x = min(posX)
        upper_left_um_y = max(posY)
        deltax, deltay = pks_from_img.get_translation()
        print(deltax,deltay)
        print(upper_left_um_x, upper_left_um_y)
        labels_res = []
        posX_res = []
        posY_res = []
        for i,x in enumerate(labels):

            x = int(deltax+(305-x_border)/2 + (posX[i] - upper_left_um_x)/0.34) # tile 2 06/09/2022;
            y = int(deltay+(154-y_border)/2  - (posY[i] - upper_left_um_y)/0.34)  # 305 and 154 are the size of tfd 512 by 256 image

            if (x > 0 ) and (y > 0):
                if (y+y_border< 2944) and (x+x_border <2866) :
                    labels_res.append (labels[i])
                    posX_res.append (x)
                    posY_res.append(y)
        print (labels_res)
        print('Number of overlapping FOVs: ',len(labels_res))
        return (labels_res, posX_res, posY_res)
    else:
        return (labels, posX, posY)
