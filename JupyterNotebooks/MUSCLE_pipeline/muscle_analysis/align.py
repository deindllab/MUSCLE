import os
import copy
import os.path 
import cv2
import tkinter as tk
import tkinter.filedialog as fd
from PIL import Image, ImageEnhance
from skimage import transform, img_as_int, exposure, img_as_ubyte
import skimage.io as io
import skimage as si
import tkinter.messagebox as mb
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import blob_log, blob_doh, blob_dog
from scipy import ndimage, spatial
from cmath import inf
import sys
import numpy as np
from scipy.optimize import curve_fit
import json
from Bio import SeqIO
from skimage.filters import gaussian
import pandas as pd
from skimage.feature import peak_local_max
from sklearn.metrics import mean_squared_error
from numpy.linalg import norm
from math import floor, ceil, exp
from skimage.registration import phase_cross_correlation
from scipy import ndimage as ndi
from scipy.io import savemat
import itertools
import scipy.spatial
from skimage.restoration import rolling_ball
from . import muscle_sequencing, pks_from_img, Polywarp, star, align, load_data
root = tk.Tk()
root.attributes("-topmost", True)
root.withdraw()

def alignment (path_smFRET, POS, ALEX, tr_R2G, apriori_tr, x_coord, y_coord, sequence, read2 ):
    """
    The function is used for the alingment between single molecules from smFRET imaging and clusters from FASTQ file.

    Args:
        path_smFRET (string): The absolute path of the folder with smFRET data.
        POS: the information about positions. 
        ALEX (boolean): If imaging was in ALEX regime.
        tr_R2G : forward transformation 
        apriori_tr: original transformation  
        x_coord: list of the x coordinates from the FASTQ file
        y_coord: list of the y coordinates from the FASTQ file
        sequence: list of the sequences from the FASTQ file, in case of paired-end sequencing, it is 2D array
        read2: in case of paired-end sequencing, this flag should be True
    Returns:
         It returns the following parameteres which have been chosen here: 
         #r - Half-width of the molecule aperture for trace extraction, i.e. for r = 3 it is -3:3
         #rb_rad - radius for the ralling ball background estimation
         And creates .mat file with "good positions"
    
    """ 
    
    if read2:
        sequence_1 = sequence[0]
        sequence_2 = sequence[1] ##for paired end
    else:
        sequence_1 = sequence
   
    print (len(sequence_1))
    print (len(sequence_2))
    if ALEX:
        green_frames = np.add(1,np.multiply(2,range(10)))
        red_frames = np.subtract(green_frames,1)
    x_border = 500
    y_border = 300
    counter = 0  
    good_pos = []
    matched_sequences = []
    matched_centers_red = []
    matched_centers_green = []
    matched_centers_FQ = []
    matched_dist = []
    processed_pos = []
    FQ_shifts = []
    FQ_shifts_pos = []
    matched_sequences_2 = []
    
    rb_rad = int(input("Please enter the radius for the ralling ball background estimation"))
    rb = int(rb_rad/2) # Half-width of the background aperture for trace extraction
    r = int(input("Please enter the half-width of the molecule aperture for trace extraction, , i.e. for r = 3 it is -3:3"))
    
    
    labels_res, posX_res, posY_res = load_data.extract_pos_info(POS, res = True)
    current_direct = fd.askdirectory(title = "Choose the output folder")
    os.makedirs(os.path.join(current_direct, 'QC_composites'), exist_ok=True)
    os.makedirs(os.path.join(current_direct, 'QC_composites_raw'), exist_ok=True)
    
    
    #r = 3 # Half-width of the molecule aperture for trace extraction, i.e. for r = 3 it is -3:3
    #rb = 5 # Half-width of the background aperture for trace extraction 
    for pos in labels_res :
        if os.path.exists(os.path.join(path_smFRET,pos,'B_Green','img_000000000.tiff')):
            path_smFRET_file = os.path.join(path_smFRET,pos,'B_Green','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_hairpinGreen','img_000000000.tiff')):
            path_smFRET_file = os.path.join(path_smFRET,pos,'B_hairpinGreen','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_Green_hairpin','img_000000000.tiff')):
            path_smFRET_file = os.path.join(path_smFRET,pos,'B_Green_hairpin','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_Green_Normal','img_000000000.tiff')):
            path_smFRET_file = os.path.join(path_smFRET,pos,'B_Green_Normal','img_000000000.tiff')
        else:
            continue

        if (not os.path.exists(os.path.join(current_direct, pos + "_traces.mat")))&(not pos in processed_pos):

            processed_pos.append(pos)
            pos_direct = os.path.join(current_direct, pos)
            os.makedirs(pos_direct, exist_ok=True)
            log_file = open(os.path.join(pos_direct, pos + "_log.txt"), 'w')
            log_file.write('Working on: '+ pos+'\n')
            print('Working on: '+ pos+'\n')
            FRET_coord = []
            seq_coord = []
            point = []
            counter = 0
            
            idx = [i for i,x in enumerate(labels_res) if labels_res[i]==pos]
            X_c = [posX_res[i] for i in idx]
            Y_c = [posY_res[i] for i in idx]

            #X_c, Y_c = muscle_sequencing.scaling_seq(X_c[0], Y_c[0]) #Misha: should be imported from muscle_sequencing.py

            point.append(X_c[0])
            point.append(Y_c[0])
            point = np.array(point)
            idx = [i for i,x in enumerate(x_coord) if x_coord[i] >= X_c[0] and (x_coord[i] <= X_c[0]+x_border) and (y_coord[i] >= Y_c[0]) and (y_coord[i] <= Y_c[0]+y_border)]
            x_FQ = [x_coord[i] for i in idx]
            y_FQ = [y_coord[i] for i in idx]
            seq_t = [sequence_1[i] for i in idx]
            if read2:
                seq_t_2 = [sequence_2[i] for i in idx]
            if len(seq_t)<100:
                log_file.write('Not enough (<100) clusters: '+ str(len(seq_t))+ '\n')
                continue


            x_FQ = np.subtract(x_FQ,X_c[0])
            y_FQ = np.subtract(y_FQ,Y_c[0])
            fastq_image1 = muscle_sequencing.generate_img(y_FQ,x_FQ, 0, 0, y_border, x_border, 1, True, 1)  #Misha: should be imported from muscle_sequencing.py
            fastq_image1.save(os.path.join(pos_direct,pos+'_FQ.png'))

            fig, ax = plt.subplots()
            ax.imshow(fastq_image1)
            plt.show()


            img_smFRET = io.imread(path_smFRET_file)

            # Averaging the first 10 frames to select peaks
            ALEX = True
            if ALEX:
#             img_t = np.mean(img_smFRET[green_frames,::], axis = 0)
                img_t = np.mean(img_smFRET[red_frames,::], axis = 0)
            else:
                img_t = np.mean(img_smFRET[0:10,::], axis = 0)
            img_t = img_t.astype("ushort")
            

            red = img_t[256:,:]
            green = img_t[:256,:]
            red = red - rolling_ball(red, radius=rb_rad)
            green = green - rolling_ball(green, radius=rb_rad)
            green = transform.warp(green,tr_R2G, preserve_range = True)
            combined = red + green # Consider adding the red excitation channel, though there are some difficulties, e.g. beads and int scaling
            fig, ax = plt.subplots()
            ax.imshow(combined)
            blobs_log = blob_log(combined, max_sigma=10, num_sigma=10, threshold=30) # Was 1000 for 19/07/2022
    #             Was 300 for 06/09/2022
            CM = []
            #r = 3
            [h,w] = red.shape
            n_frames = img_smFRET.shape[0]
            FQ_centers = np.concatenate((x_FQ, y_FQ)).reshape((-1, 2), order='F')

            for i, blob in enumerate(blobs_log):
                x, y, d = blob
                if x>r and x<(h-r) and y>r and y<(w-r):
                    temp = ndimage.measurements.center_of_mass(combined[int(x-r):int(x+r+1),int(y-r):int(y+r+1)])
                    CM.append(np.flip(np.add(temp, [x-r,y-r])))

                    c = plt.Circle(CM[-1], 3, color="red", linewidth=1, fill=False)
                    ax.add_patch(c)
            ax.set_axis_off()
            plt.savefig(os.path.join(pos_direct, pos + "_smFRET_peaks.tif"))
            plt.show()
            smFRET_centers = np.array(CM)
    #         res, idx = count_nearest_pts(rough_tr(movie_centers), seq_centers, 8)
    #         movie_centers1 = movie_centers[idx[np.where(res != inf)]]
    #         seq_centers1 = seq_centers[np.where(res != inf)]

    #         if not beadsOK:
                # Calculating the shift between the FASTQ and cluster images using cross-correlation
            smFRET_centers_tr = apriori_tr(smFRET_centers)
            smFRET_centers_tr_image_array = np.asarray(muscle_sequencing.generate_img(smFRET_centers_tr[:,1],smFRET_centers_tr[:,0], 0, 0, y_border, x_border, 1, True, 1))
            fastq_image_lib_array = np.asarray(muscle_sequencing.generate_img(y_FQ,x_FQ, 0, 0, y_border, x_border, 1, True, 1))
            fig, ax = plt.subplots()
            ax.imshow(smFRET_centers_tr_image_array)
            plt.show()
            image_smFRET_tr = Image.fromarray(smFRET_centers_tr_image_array)
            image_smFRET_tr.save(os.path.join(pos_direct, pos + "_smFRET_peaks_tr.tif"))

            fig, ax = plt.subplots()
            ax.imshow(fastq_image_lib_array)
            plt.show()
            image_smFRET_tr = Image.fromarray(fastq_image_lib_array)
            image_smFRET_tr.save(os.path.join(pos_direct, pos + "_FQ_peaks.tif"))

            shift = phase_cross_correlation(fastq_image_lib_array,smFRET_centers_tr_image_array)
            shift = shift[0]
            print(shift)

            rough_tr1 = copy.deepcopy(apriori_tr)
            dx = apriori_tr.params[0,2] + shift[1] # Check sign!
            dy = apriori_tr.params[1,2] + shift[0]
            rough_tr1.params[0,2] = dx
            rough_tr1.params[1,2] = dy

            np.save(os.path.join(pos_direct, pos + "_auto_bead_tr"), rough_tr1)

            smFRET_centers_tr = rough_tr1(smFRET_centers)
            smFRET_centers_tr_image_array = np.asarray(muscle_sequencing.generate_img(smFRET_centers_tr[:,1],smFRET_centers_tr[:,0], 0, 0, y_border, x_border, 1, True, 1))
            fastq_image_lib_array = np.asarray(muscle_sequencing.generate_img(y_FQ,x_FQ, 0, 0, y_border, x_border, 1, True, 1))
            fig, ax = plt.subplots()
            ax.imshow(smFRET_centers_tr_image_array)
            plt.show()
            image_smFRET_tr = Image.fromarray(smFRET_centers_tr_image_array)
            image_smFRET_tr.save(os.path.join(pos_direct, pos + "_smFRET_peaks_upd_tr.tif"))        

            res, idx = pks_from_img.count_nearest_pts(rough_tr1(smFRET_centers), FQ_centers, 2)
            movie_centers1 = smFRET_centers[idx[np.where(res != inf)]]
            seq_centers1 = FQ_centers[np.where(res != inf)]

            if seq_centers1.shape[0]<4:
                print('Not enough smFRET peaks (<4)')
                log_file.write('Not enough smFRET peaks (<4)')
                continue

            tr = transform.estimate_transform("similarity", src=movie_centers1, dst=seq_centers1)
            tr_inv = transform.estimate_transform("similarity", src=seq_centers1, dst=movie_centers1)
            np.save(os.path.join(pos_direct, pos + "_final_tr"), tr)
            np.save(os.path.join(pos_direct, pos + "_final_tr_inv"), tr_inv)
            log_file.write('Matched clusters before transform update: ' + str(seq_centers1.shape[0])+'\n')
            print('Matched clusters before transform update: ',seq_centers1.shape[0])

            res, idx = pks_from_img.count_nearest_pts(tr(smFRET_centers), FQ_centers, 3) # Was 4 for 19/07/2022
            centers_matched = smFRET_centers[idx[np.where(res != inf)]]
            idx_t = np.where(res != inf)
            seq_matched = [seq_t[i] for i in idx_t[0]]
            if read2:
                seq_matched_2 = [seq_t_2[i] for i in idx_t[0]]
            FQ_centers_matched = FQ_centers[idx_t]
            distances_matched = res[idx_t]

            centers_red = centers_matched
            centers_green = tr_R2G(centers_matched)
            #r = 3
            # Weed out positions that are too close to the edge
            idx_t = np.where((centers_red[:,0]>rb) & (centers_red[:,0]<(w-rb)) & (centers_red[:,1]>rb) & (centers_red[:,1]<(h-rb))  & 
                             (centers_green[:,0]>rb) & (centers_green[:,0]<(w-rb)) & (centers_green[:,1]>rb) & (centers_green[:,1]<(h-rb)))
            centers_red = centers_red [idx_t]
            centers_green = centers_green [idx_t]
            distances_matched = distances_matched[idx_t]
            seq_matched = [seq_matched[i] for i in idx_t[0]]
            if read2:
                seq_matched_2 = [seq_matched_2[i] for i in idx_t[0]]
            FQ_centers_matched = FQ_centers_matched[idx_t]

            matched_sequences.append(seq_matched)
            if read2:
                matched_sequences_2.append(seq_matched_2)
            matched_centers_red.append(centers_red)
            matched_centers_green.append(centers_green)
            matched_dist.append(distances_matched)
            matched_centers_FQ.append(FQ_centers_matched)
            good_pos.append(pos)
            log_file.write('Matched clusters after transform update: ' + str(len(seq_matched))+'\n')
            log_file.write('Out of ' + str(len(smFRET_centers)) + ' smFRET peaks and ' + str(len(FQ_centers)) + ' clusters'+'\n')
            log_file.write('Percentage of matched smFRET peaks: ' + str(int(100*len(seq_matched)/len(smFRET_centers)))+'\n')
            log_file.write('Percentage of matched clusters: ' + str(int(100*len(seq_matched)/len(FQ_centers)))+'\n')
            print('Matched clusters after transform update: ',len(seq_matched))
            print('Out of ',len(smFRET_centers), ' smFRET peaks and ', len(FQ_centers), ' clusters')
            print('Percentage of matched smFRET peaks: ',int(100*len(seq_matched)/len(smFRET_centers)))
            print('Percentage of matched clusters: ',int(100*len(seq_matched)/len(FQ_centers)))
            # Saving transformed bead and molecule images for QC
    #         img_beads_tfd = Image.fromarray(transform.warp(img_beads_F, tr_inv, output_shape = [y_border, x_border]))
    #         img_beads_tfd.save(os.path.join(pos_direct, pos + "_beads_transformed.tif"))
            img_smFRET_tfd = Image.fromarray(transform.warp(combined, tr_inv, output_shape = [y_border, x_border]))
            img_smFRET_tfd.save(os.path.join(pos_direct, pos + "_smFRET_transformed.tif"))
            centers_red_tfd = tr(centers_red)
            img_centers_red_tfd = muscle_sequencing.generate_img(centers_red_tfd[:,1],centers_red_tfd[:,0], 0, 0, y_border, x_border, 1, True, 1)
            img_centers_red_tfd.save(os.path.join(pos_direct, pos + "_smFRET_peaks_matched_transformed.tif"))

            centers_red_tfd = tr(smFRET_centers)
            img_centers_red_tfd = muscle_sequencing.generate_img(centers_red_tfd[:,1],centers_red_tfd[:,0], 0, 0, y_border, x_border, 1, True, 1)
            img_centers_red_tfd.save(os.path.join(pos_direct, pos + "_smFRET_peaks_all_transformed.tif"))

            img_centers_red_tfd = muscle_sequencing.generate_img(FQ_centers_matched[:,1],FQ_centers_matched[:,0], 0, 0, y_border, x_border, 1, True, 1)
            img_centers_red_tfd.save(os.path.join(pos_direct, pos + "_FQ_lib_matched.tif"))

            img_centers_red_tfd = muscle_sequencing.generate_img(FQ_centers[:,1],FQ_centers[:,0], 0, 0, y_border, x_border, 1, True, 1)
            img_centers_red_tfd.save(os.path.join(pos_direct, pos + "_FQ_lib_all.tif"))

    #         Saving a composite image for quick QC
            fig, ax = plt.subplots()
            QC_composite_array = np.zeros([y_border,x_border,3], dtype=np.uint8)
            red_channel = np.asarray(muscle_sequencing.generate_img(centers_red_tfd[:,1],centers_red_tfd[:,0], 0, 0, y_border, x_border, 1, True, 1))
            percentiles = np.percentile(red_channel[75:125,:], (0.5, 99.5))        
            red_channel = exposure.rescale_intensity(red_channel, in_range=tuple(percentiles))
            red_channel = img_as_ubyte(red_channel)

            green_channel = np.asarray(muscle_sequencing.generate_img(FQ_centers[:,1],FQ_centers[:,0], 0, 0, y_border, x_border, 1, True, 1))
            percentiles = np.percentile(green_channel, (0.5, 99.5))
            green_channel = exposure.rescale_intensity(green_channel, in_range=tuple(percentiles))
            green_channel = img_as_ubyte(green_channel)

            blue_channel = np.asarray(muscle_sequencing.generate_img(FQ_centers_matched[:,1],FQ_centers_matched[:,0], 0, 0, y_border, x_border, 1, True, 1))
            percentiles = np.percentile(blue_channel, (0.5, 99.5))
            blue_channel = exposure.rescale_intensity(blue_channel, in_range=tuple(percentiles))
            blue_channel = img_as_ubyte(blue_channel)

            QC_composite_array[:,:,0] = red_channel
            QC_composite_array[:,:,1] = green_channel
            QC_composite_array[:,:,2] = blue_channel

            QC_composite = Image.fromarray(QC_composite_array, mode="RGB")
            QC_composite.save(os.path.join(current_direct, 'QC_composites_raw' , pos + "_QC.tif"))

            ax.imshow(QC_composite)
            for blob in FQ_centers_matched:
                c = plt.Circle(blob, 3, color="red", linewidth=1, fill=False)
                ax.add_patch(c)
            ax.set_axis_off()
            plt.savefig(os.path.join(current_direct, 'QC_composites' , pos + "_QC.png"))
            plt.show()


    # Save good pos data to be able to extract traces later
    
    #frame_rate = 5
    mdict1 = {
        "good_pos": good_pos,
        "matched_sequences": matched_sequences,
        "matched_sequences_2": matched_sequences_2,
        "matched_centers_red": matched_centers_red,
        "matched_centers_green": matched_centers_green,
        "matched_distances": matched_dist,
        "matched_FASTQ": matched_centers_FQ,
        "FASTQ_shifts": FQ_shifts
    }

    savemat(os.path.join(current_direct, "good_pos.mat"), mdict1)
    
    print('Number of matched FOVs: ',len(good_pos))
    return (rb_rad, r)