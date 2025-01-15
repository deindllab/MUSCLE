# constants.py
from skimage import transform
import numpy as np

# Define the positions of the red and green channel images within the camera frame
# See Define_channel.jpg for more details

# r_x_start: starting image coordinate of red channel x axis
# r_x_end: ending image coordinate of red channel x axis
# r_y_start: starting image coordinate of red channel y axis
# r_y_end: ending image coordinate of red channel y axis
# g_x_start: starting image coordinate of green channel x axis
# g_x_end: ending image coordinate of green channel x axis
# g_y_start: starting image coordinate of green channel y axis
# g_y_end: ending image coordinate of green channel y axis

# NOTE!: sizes of the red and green channels are expected to be the same!

red_x_start = 0
red_x_end = 512
red_y_start = 256
red_y_end = 512
green_x_start = 0
green_x_end = 512
green_y_start = 0
green_y_end = 256

# Defining the relationship between the microscope stage coordinate system (from the position list)
# and the microscope image coordinate system. We are assuming that only 90-degree rotations and 
# reflections are required to align the two coordinate systems. See stage_cor.jpg for more details
        
stage_rotation = 1 # determine whether the stage axes x and y are inverted compare to a standard coordinate system. -1 for inverted axis, 1 for normal axis.
stage_x_cor = 1 # -1 if the FOV images are taken from minimum value to 0 in position file, 1 for 0 to maximum value
stage_y_cor = -1 # 1 if the FOV images are taken from minimum value to 0 in position file, -1 for 0 to maximum value

# Defining the relationship between the FASTQ coordinate system (from the position list)
# and the microscope image coordinate system. We are assuming that only 90-degree rotations and 
# reflections are required to align the two coordinate systems.
        
FQ_rotation = -1 # determine whether the FASTQ axes x and y are inverted compare to a standard coordinate system. -1 for inverted axis, 1 for normal axis.
FQ_x_cor = -1 # sign of the FASTQ x axis direction vs that of the smFRET image
FQ_y_cor = -1 # sign of the FASTQ x axis direction vs that of the smFRET image

# Define pixel sizes of the single-molecule FRET microscope and the MiSeq sequencer in nm/pixel

smFRET_pixel = 204 # nm/pixel
MiSeq_pixel = 340 # nm/pixel

tile_size_y = 2944 # pixel
tile_size_x = 2866 # pixel

# Canvas size for the FASTQ image
max_y = 3000 # pixel
max_x = 3000 # pixel

# The size of the 
template_size = 800 # pixel

converted_x = (red_x_end - red_x_start)*smFRET_pixel/MiSeq_pixel
converted_y = (red_y_end - red_y_start)*smFRET_pixel/MiSeq_pixel
x_border = int(1.5*converted_x)
y_border = int(1.5*converted_y)


# Defining the apriori transformations
scale = smFRET_pixel/ MiSeq_pixel
apriori_tr_original = transform.SimilarityTransform(scale = scale)
apriori_tr = transform.SimilarityTransform(scale = scale)

# Centering the tranformation to the x_border and y_border
delta = 0.5*np.subtract([x_border,y_border],apriori_tr([red_x_end-red_x_start,red_y_end-red_y_start]))
dx, dy = delta[0]
apriori_tr.params[0,2] = dx
apriori_tr.params[1,2] = dy
apriori_tr_inv = transform.SimilarityTransform(scale = 1/scale)

# Centering the tranformation to the x_border and y_border
delta = 0.5*np.subtract([red_x_end-red_x_start,red_y_end-red_y_start],apriori_tr_inv([x_border,y_border]))
dx, dy = delta[0]
apriori_tr_inv.params[0,2] = dx
apriori_tr_inv.params[1,2] = dy

# Green and red frames for ALEX
green_frames = np.add(1,np.multiply(2,range(10)))
red_frames = np.subtract(green_frames,1)