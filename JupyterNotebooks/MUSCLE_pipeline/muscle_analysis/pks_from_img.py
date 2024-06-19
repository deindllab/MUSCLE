import numpy as np
import os
import matplotlib.pyplot as plt
from skimage.feature import blob_log, blob_doh, blob_dog
from scipy import ndimage, spatial
from skimage.restoration import rolling_ball
import skimage.io as io
from skimage import transform, img_as_int, exposure, img_as_ubyte
from . import muscle_sequencing, load_data
from PIL import Image
import tkinter.filedialog as fd
from cmath import inf

def blob_detection(img, min_sigma, max_sigma, threshold, method=0):
    """This function is mostly used for detecting the beads in any image.

    Args:
        img_path (string): The absolute path of the input image.
        min_sigma (int): The minimum sigma, lower it is, smaller the blob will be detected.
        max_sigma (int): The maximum sigma, higher it is, bigger the blob will be detected.
        threshold (float): Higher it is, higher the intensities of blobs.
        method (int, optional): 0 for Difference of Gaussian (DoG) and 1 for Determinant of Hessian (DoH). 
        They should be applied with different combination of parameters. DoG is more suitable for fret movies,
        while DoH is more suitable for sequencing images. Defaults to 0.

    Returns:
        centers: A numpy array containing the coordinates of all the centers.
    """
    #img = io.imread(img_path)
    fig, ax = plt.subplots()
    ax.imshow(img)
    if method == 0:
        blob = blob_dog(
            img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold
        )
    else:
        blob = blob_doh(
            img, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold
        )
    i = 0
    r = 3
    centers = []
    h, w = img.shape
    for blob in blob:
        y, x, r = blob
        if y > r and y < (h - r) and x > r and x < (w - r):
            centers.append(
                ndimage.measurements.center_of_mass(
                    img[int(y - r) : int(y + r + 1), int(x - r) : int(x + r + 1)]
                )
            )
            centers[i] = list(np.add(np.flip(centers[i]), [x - r, y - r]))
            x1, y1 = centers[i]
            c = plt.Circle([x1, y1], 3, color="red", linewidth=1, fill=False)
            ax.add_patch(c)
            
            i += 1
    ax.set_axis_off()
    plt.show()
    return np.array(centers)

def count_nearest_pts(src, dst, radius):
    """Counting the number of nearest neighbors for each given point.

    Args:
        src (numpy array): (N, 2) shape array. Build the kd tree based on this.
        dst (numpy array): (N, 2) shape array. For each point in this array, find the nearest neighbors in src array.
        radius (int): The maximum searching radius.

    Returns:
        res, idx: res is the distance for the point and its neighbor, 'inf' means no neighbor in given search radius. 
        idx is the index for the neighbor in src array.
    """
    tree = spatial.KDTree(src)
    res, idx = tree.query(dst, k=1, distance_upper_bound=radius)
    for i in range(0, len(idx)):
        idx_t = np.argwhere(idx == idx[i])
        if len(idx_t) > 1:
            res_t = [res[j] for j in idx_t]
            if res[i] > min(res_t): res[i] = inf
    return res, idx

def combined_image(path_smFRET, POS, ALEX, tr_R2G, apriori_tr_original):
    """Creating the combined image from the sm_FRET FOVs

    Args:
        path_smFRET: path to the folder with smFRET movies.
        POS: position list.
        ALEX (bool): was the data collected using alternating laser excitation (ALEX)?.
        tr_R2G: the polynomial reverse transformation for mapping the Cy3 and Cy5 channels
        apriori_tr_original: original apriori transformation from smFRET to sequencing coordinates
    Returns:
        The function creates an image that combines all smFRET FOVs in order to perform preliminary alignment with FASTQ tiles
    """
    size1 = apriori_tr_original([[256,512]]).astype(int)
    labels, posX, posY  = load_data.extract_pos_info(POS) #from module load_data_extracting the coordinates and the label from the position, if False, then just extract, if True choose matching 
    minx = min(posX)
    maxx = max(posX)
    miny = min(posY)
    maxy = max(posY)
    size = [maxy-miny+50+size1[0,0], maxx-minx+50+size1[0,1]]
    size = np.divide(size,0.34).astype(int)
    img = np.zeros(size)
    counter = 1
    
    
    show_im = True
    if ALEX:
        green_frames = np.add(1,np.multiply(2,range(10)))
        red_frames = np.subtract(green_frames,1)
    for pos in labels:
        if os.path.exists(os.path.join(path_smFRET,pos,'B_Green','img_000000000.tiff')):
            path1 = os.path.join(path_smFRET,pos,'B_Green','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_hairpinGreen','img_000000000.tiff')):
            path1 = os.path.join(path_smFRET,pos,'B_hairpinGreen','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_Green_hairpin','img_000000000.tiff')):
            path1 = os.path.join(path_smFRET,pos,'B_Green_hairpin','img_000000000.tiff')
        elif os.path.exists(os.path.join(path_smFRET,pos,'B_Green_Normal','img_000000000.tiff')):
            path1 = os.path.join(path_smFRET,pos,'B_Green_Normal','img_000000000.tiff')
        else:
            continue


        print(pos)
        img_smFRET = io.imread(path1)

        # Averaging the first 10 frames to select peaks

        if ALEX:
#         img_t = np.mean(img_smFRET[green_frames,::], axis = 0)
            img_t = np.mean(img_smFRET[red_frames,::], axis = 0)
        else:
            img_t = np.mean(img_smFRET[0:10,::], axis = 0)
        img_t = img_t.astype("ushort")
        rb_rad = 10


        red = img_t[256:,:]
        green = img_t[:256,:]
        red = red - rolling_ball(red, radius=rb_rad)
        green = green - rolling_ball(green, radius=rb_rad)
        green = transform.warp(green,tr_R2G, preserve_range = True)
        combined = red + green # Consider adding the red excitation channel, though there are some difficulties, e.g. beads and int scaling
        
        if show_im:
            fig, ax = plt.subplots()
            ax.imshow(combined)
        blobs_log = blob_log(combined, max_sigma=10, num_sigma=10, threshold=30) # Was 1000 for 19/07/2022
    #             Was 300 for 06/09/2022
        CM = []
        r = 3
        [h,w] = red.shape

        for i, blob in enumerate(blobs_log):
            x, y, d = blob
            if x>r and x<(h-r) and y>r and y<(w-r):
                temp = ndimage.measurements.center_of_mass(combined[int(x-r):int(x+r+1),int(y-r):int(y+r+1)])
                CM.append(np.flip(np.add(temp, [x-r,y-r])))
                if show_im:
                    c = plt.Circle(CM[-1], 3, color="red", linewidth=1, fill=False)
                    ax.add_patch(c)
        if show_im:
            ax.set_axis_off()
    #     plt.savefig(os.path.join(pos_direct, pos + "_smFRET_peaks.tif"))
            plt.show()
        smFRET_centers = np.array(CM)
        smFRET_centers_tr = apriori_tr_original(smFRET_centers)
        smFRET_centers_tr_image_array = np.asarray(muscle_sequencing.generate_img(smFRET_centers_tr[:,1],smFRET_centers_tr[:,0], 0, 0,size1[0,0], size1[0,1], 1, True, 1))
   
        x1 = [posX[i] for i,x in enumerate(labels) if labels[i]==pos]
        y1 = [posY[i] for i,x in enumerate(labels) if labels[i]==pos]
        id1 = int((maxy-y1[0])/0.34)
        id2 = int((x1[0]-minx)/0.34)

        img[id1:id1+size1[0,0],id2:id2+size1[0,1]] = smFRET_centers_tr_image_array
        counter = counter + 1
    img2 = Image.fromarray(img.astype(np.byte))
    img2 = img2.convert('L')
    current_direct = fd.askdirectory(title = "Choose the output folder")
    img2.save(current_direct+'/Combined.png')
    
    print("combined image is successfully saved")
    
    
    
    
    
    
    
    
    
   

    

