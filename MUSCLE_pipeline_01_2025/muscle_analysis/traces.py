import numpy as np 
import os 
import os.path
import skimage.io as io
from math import floor, ceil, exp
import tkinter.filedialog as fd
import skimage as si
from scipy.io import savemat, loadmat
from .constants import *

def extract_traces(path_smFRET,read2, rb_rad, r):
    """
    This function is creating the .mat files containing matched smFRET traces and sequences for molecules from FOVs that were aligned in the previous step

    Args:
        path_smFRET - path to the folder containing the smFRET movies. 
        read2 - flag for paired-end sequencing
    Returns:
        rb_rad - rolling ball radius
        r - Half-width of the molecule aperture for trace extraction, i.e. for r = 3 it is -3:3
      
        "time": time information from smFRET trace 
        "Cy3": Cy3 intensity 
        "Cy5": Cy5 intensity 
        "Seq": sequence based on the read1,
        "Seq_2": sequence based on read2 (if not paired-end, it's an empty array),
        "Dist": distance betweeen a cluster and a corresponding single molecule,
        "x": x coordinates of single molecules from the Cy5 channel,
        "y": y coordinates of single molecules from the Cy5 channel,
        "x_FQ": x coordinates of respective clusters from the fastq file,
        "y_FQ": y coordinates of respective clusters from the fastq file
             
    """
    
    file_path = fd.askopenfilename(title = "Choose the good_pos.mat file") 
    mdict1 = loadmat(file_path)
    
    current_direct = fd.askdirectory(title = "Choose the output folder")


    # Extracting traces, median BG version
    frame_rate = int(input("Please, choose the frame rate(Hz) "))
    rb = int(rb_rad/2)  # Half-width of the background aperture for trace extraction 
    # (bacground is calculated as a median of intensities in a rectangular aperture between r and rb)
    # a = np.subtract(range(2*rb+1), rb)
    
    at = [[j,i] for i in range(2*rb+1) for j in range(rb-r)]
    bt = [[j,i] for i in range(2*rb+1) for j in range(r+rb+1,2*rb+1)]
    ct = [[j,i] for i in range(rb-r) for j in range(rb-r,r+rb+1)]
    dt = [[j,i] for i in range(r+rb+1,2*rb+1) for j in range(rb-r,r+rb+1)]
    idx_bg = np.subtract(np.concatenate([at,bt,ct,dt]), [rb,rb]) # Indexes that define the background aperture

    for i, pos in enumerate(mdict1["good_pos"]):
        print(pos)
        path_smFRET_file = []
        for current_dir, dirs, files in os.walk(os.path.join(path_smFRET,pos)) :
            for el in files:
                if (el.split('.')[-1].lower() == 'tif') or (el.split('.')[-1].lower() == 'tiff'):
                    path_smFRET_file = os.path.join(current_dir,el)
                    break
            else: continue
            break
        if path_smFRET_file:
            img_smFRET = io.imread(path_smFRET_file)
            print('Working on: '+pos)
            
            seq_matched = mdict1 ["matched_sequences"][0][i]
            if read2:
                seq_matched_2 =  mdict1 ["matched_sequences_2"][0][i]
            else:
                 seq_matched_2 = []
            dist_matched =  mdict1 ["matched_distances"][0][i]
            centers_red = mdict1 ["matched_centers_red"][0][i]
            centers_green = mdict1 ["matched_centers_green"][0][i]
            centers_FQ = mdict1["matched_FASTQ"][0][i]
            n_traces = len(centers_red)
            [h,w] = [red_y_end-red_y_start,red_x_end-red_x_start]
            n_frames = img_smFRET.shape[0]
            print (i)
            weights_red = np.zeros([n_traces,2*r+1,2*r+1])
            weights_green = np.zeros([n_traces,2*r+1,2*r+1])
            traces_red =  np.zeros([n_traces,n_frames])
            traces_green = np.zeros([n_traces,n_frames])
            # Extracting the traces from matched smFRET peaks
            # Extracting the traces from matched smFRET peaks
            # Select matching peaks from peak_locations and their sequences


            for j,coord in enumerate(centers_red):
                x, y = coord
                x0 = int(x-r)
                x1 = int(x+r+1)
                y0 = int(y-r)
                y1 = int(y+r+1)
                dx = x-x0
                dy = y-y0

                it = np.nditer(weights_red[j], flags=['multi_index'], op_flags=['readwrite'])
                for w1 in it:
                    yt,xt = it.multi_index
                    w1[...] = 2*exp(-0.4*((xt-dx)**2+(yt-dy)**2))
                it.close()

                x, y = centers_green[j]
                x0 = int(x-r)
                x1 = int(x+r+1)
                y0 = int(y-r)
                y1 = int(y+r+1)
                dx = x-x0
                dy = y-y0

                it = np.nditer(weights_green[j], flags=['multi_index'], op_flags=['readwrite'])
                for w1 in it:
                    yt,xt = it.multi_index
                    w1[...] = 2*exp(-0.4*((xt-dx)**2+(yt-dy)**2))
                it.close()

            for k in range(n_frames):
                if k%100 == 0: print("Working on frame "+str(k))
                red1 = img_smFRET[k,red_y_start:red_y_end,red_x_start:red_x_end]
                red1 = red1 - si.restoration.rolling_ball(red1, radius=rb_rad)
                green1 = img_smFRET[k,green_y_start:green_y_end,green_x_start:green_x_end]
                green1 = green1 - si.restoration.rolling_ball(green1, radius=rb_rad)

                for j,coord in enumerate(centers_red):
                    x, y = coord
                    x0 = int(x-r)
                    x1 = int(x+r+1)
                    y0 = int(y-r)
                    y1 = int(y+r+1)
                    dx = x-x0
                    dy = y-y0


                    index = np.add(idx_bg,[int(x), int(y)])
              
                    BG = np.median(red1[index[:,1], index[:,0]])
                    red2 = np.subtract(red1[y0:y1,x0:x1], BG)
                    traces_red[j,k] = np.sum(np.multiply(weights_red[j],red2))


                    x, y = centers_green[j]
                    x0 = int(x-r)
                    x1 = int(x+r+1)
                    y0 = int(y-r)
                    y1 = int(y+r+1)
                    dx = x-x0
                    dy = y-y0

                    index = np.add(idx_bg,[int(x), int(y)])
        
                    BG = np.median(green1[index[:,1], index[:,0]])
                    green2 = np.subtract(green1[y0:y1,x0:x1], BG)
                    traces_green[j,k] = np.sum(np.multiply(weights_green[j],green2))

           
            # Save traces and sequences
            
            mdict = {
            "time": np.divide(range(n_frames),frame_rate),
            "Cy3": traces_green,
            "Cy5": traces_red,
            "Seq": seq_matched,
            "Seq_2": seq_matched_2,
            "Dist": dist_matched,
            "x": centers_red[:,0],
            "y": centers_red[:,1],
            "x_FQ": centers_FQ[:,0],
            "y_FQ": centers_FQ[:,1]
               
            }
            
            savemat(os.path.join(current_direct, pos + "_traces.mat"), mdict)
            print("traces saved")
