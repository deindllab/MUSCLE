# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:48:42 2022

@author: Anton
"""
import warnings

import tkinter as tk
from tkinter import filedialog as fd
import os
import os.path 
from skimage import transform
import skimage.io as io
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import blob_log
# from scipy import ndimage

from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt
# import skimage as si
# import skimage.io as io
# from skimage.feature import blob_log
# from scipy import ndimage
# import matplotlib.pyplot as plt
# from matplotlib import axes
# from matplotlib import figure
# import lmfit
from numpy import exp, power
# from lmfit import Minimizer, Parameters
from numba import jit
from skimage.restoration import rolling_ball
import concurrent.futures

@jit
def _gaussian(M, A, sigma, bg, x0, y0):
    x, y = M
    return bg + A*exp(-0.5*power(sigma,-2)*(power(x-x0,2) + power(y-y0,2)))

@jit
def _jacobian(M, A, sigma, bg, x0, y0):
    x, y = M
    e = exp(-0.5*power(sigma,-2)*(power(x-x0,2) + power(y-y0,2)))
    aes =  A*e*power(sigma,-2)
    return np.transpose(np.array([e, A*e*power(sigma,-3)*(power(x-x0,2)+power(y-y0,2)), np.ones(len(x)),aes*(x-x0), aes*(y-y0)]))

rb_rad = 10

# @jit
def stack_subtract_background(stack, n_frames):
    for i in range(n_frames):
        stack[i,:,:] = stack[i,:,:] - rolling_ball(stack[i,:,:], radius = rb_rad) 
    return stack

# @jit
def stack_fit(stack, x0, y0, sigma, n_frames):
    r = 3
    X = range(2*r+1)
    Y = range(2*r+1)
    Y , X = np.meshgrid(Y, X)
    X = X.ravel()
    Y = Y.ravel()
    xdata = np.vstack((X, Y))
    def _jacobian_fix(M, A, bg):
        x, y = M
        return np.transpose(np.array([exp(-0.5*power(sigma,-2)*(power(x-x0,2) + power(y-y0,2))), np.ones(len(x))]))
    def _gaussian_fix(M, A, bg):
        x, y = M
        return bg + A*exp(-0.5*power(sigma,-2)*(power(x-x0,2) + power(y-y0,2)))
    
    out = np.zeros(n_frames)
    for i in range(n_frames):
        frame = stack[i,:,:]
        minim = np.min(frame)
        p0_fix = [np.max(frame)-minim, minim]
        try:
            popt = curve_fit(_gaussian_fix, xdata, frame.ravel(), p0_fix, jac = _jacobian_fix)
            out[i] = popt[0][0]
        except:
            out[i] = -10000
    return out

if __name__ == '__main__':
    warnings.warnings("ignore")
    root = tk.Tk()
    root.attributes("-topmost", True)
    root.withdraw()
    
    #Instead of calculating the bead mapping we can load the transformation generated with channel_map
    tr_G2R = transform.PolynomialTransform()
    file_path = fd.askopenfilename(title = "Choose the forward transform file", initialdir = "C:/Users/Anton/Documents/Jupyter home/ExampleData/Test_folder/220717_FC_Nano_200nm_Multicolor/original")
    tr_G2R.params = np.load(file_path)
    tr_R2G = transform.PolynomialTransform()
    file_path = fd.askopenfilename(title = "Choose the inverse transform file", initialdir = os.path.dirname(file_path))
    tr_R2G.params = np.load(file_path)
    
    file_path = fd.askopenfilename(title = "Choose the movie")
    rb_rad = 10
    r = 3
    p0 = [0, 2.0, 0, 3.0, 3.0]
    img_smFRET = io.imread(file_path)
    # Averaging the first 10 frames to select peaks
    img_t = np.mean(img_smFRET[0:10,::], axis = 0)
    img_t = img_t.astype("ushort")
    rb_rad = 10
    #print(img1[10:20,10:20])
    #img = np.zeros([20,20])
    #img[10,5] = 1000
    red = img_t[256:,:]
    green = img_t[:256,:]
    red = red - rolling_ball(red, radius=rb_rad)
    green = green - rolling_ball(green, radius=rb_rad)
    green = transform.warp(green,tr_R2G, preserve_range = True)
    combined = red + green # Consider adding the red excitation channel, though there are some difficulties, e.g. beads and int scaling
    fig, ax = plt.subplots()
    ax.imshow(combined)
    blobs_log = blob_log(combined, max_sigma=10, num_sigma=10, threshold=1000) # Was 1000 for 19/07/2022
    #             Was 300 for 06/09/2022
    CM = []
    sig = []
#             xy_init = []
    [h,w] = red.shape
    n_frames = img_smFRET.shape[0]
    X = range(2*r+1)
    Y = range(2*r+1)
    Y , X = np.meshgrid(Y, X)
    X = X.ravel()
    Y = Y.ravel()
    xdata = np.vstack((X, Y))
    for i, blob in enumerate(blobs_log):
        x, y, d = blob
        if x>r and x<(h-r) and y>r and y<(w-r):
            img1 = combined[int(x-r):int(x+r+1),int(y-r):int(y+r+1)]
#                     img1 = img1.ravel()
#                     params = model.make_params(bg = min(img1), xc = r, yc = r, sigma = 3, amplitude = max(img1)-min(img1))
#                     params["A"].value = np.max(img1)-np.min(img1)
#                     params["bg"].value = np.min(img1)
            p0[0] = np.max(img1)-np.min(img1)
            p0[2] = np.min(img1)
            try:
                popt = curve_fit(_gaussian, xdata, img1.ravel(), p0, jac = _jacobian, maxfev = 3000)
            except: continue
            popt = popt[0]
#                     min2 = Minimizer(residuals, params, fcn_args=(X,Y), fcn_kws={'z': img1.ravel()})
#                     result = min2.leastsq(Dfun=Jacobian, col_deriv=1)
#                     result = result.params                    
#                     result = result.params
            x1 = popt[3] + x - r
            y1 = popt[4] + y - r
#                     xy_init.append(popt[3], popt[4])
            CM.append([y1,x1])
            sig.append(popt[1])
            c = plt.Circle(CM[-1], 3, color="red", linewidth=1, fill=False)
            ax.add_patch(c)
    ax.set_axis_off()
#             plt.savefig(os.path.join(pos_direct, pos + "_smFRET_peaks.tiff"))
    plt.show()
    CM = np.array(CM)
    centers_red = CM
    centers_green = tr_R2G(CM)
    idx_t = np.where((centers_red[:,0]>r) & (centers_red[:,0]<(w-r)) & (centers_red[:,1]>r) & (centers_red[:,1]<(h-r))  & 
                     (centers_green[:,0]>r) & (centers_green[:,0]<(w-r)) & (centers_green[:,1]>r) & (centers_green[:,1]<(h-r)))
    centers_red = centers_red [idx_t]
    centers_green = centers_green [idx_t]
    n_traces = len(centers_red)
    [h,w] = [256,512]
    n_frames = img_smFRET.shape[0]
    traces_red =  np.zeros([n_traces,n_frames])
    traces_green = np.zeros([n_traces,n_frames])
    red1 = img_smFRET[:,256:,:]
    green1 = img_smFRET[:,:256,:]
    pool = []
    for j,coord in enumerate(centers_red):
        if j == 20: break
        print("Working on center "+str(j))
        x, y = coord
        x0 = int(x-r)
        x1 = int(x+r+1)
        y0 = int(y-r)
        y1 = int(y+r+1)
        sigma = sig[j]
        red2 = red1[:,y0:y1,x0:x1]
        y0 = x-x0
        x0 = y-y0
        sigma = sig[j]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            pool.append(executor.submit(stack_fit,stack = red2, x0 = x0, y0 = y0, sigma = sigma,n_frames =  n_frames))
        x, y = centers_green[j]
        x0 = int(x-r)   
        x1 = int(x+r+1)
        y0 = int(y-r)
        y1 = int(y+r+1)
    
        green2 = green1[:,y0:y1,x0:x1]
        y0 = x-x0
        x0 = y-y0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            pool.append(executor.submit(stack_fit, stack = green2, x0 = x0, y0 = y0, sigma = sigma,n_frames =  n_frames))

    # Save traces and sequences
    for i in concurrent.futures.as_completed(pool): print(i.result())
    
    
