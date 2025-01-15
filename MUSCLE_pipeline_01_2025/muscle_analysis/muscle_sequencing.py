import numpy as np
from skimage.filters import gaussian
from PIL import Image, ImageEnhance
from Bio import SeqIO
import tkinter.filedialog as fd
import os.path 

from .constants import *

def library_index(template, strings, min_matches):
       
    row_sums = [sum(a == b for a, b in zip(row, template)) for row in strings]
    # Create the index of elements from Seq where row_sums is above the threshold
    index = [i for i, row_sum in enumerate(row_sums) if row_sum > min_matches]
    return index

def get_pos(record):
    
    """This function is used to extract the information from the single record of the fastq file.

    Args: record - record from the .fastq file

    Returns:
        tile_num - number of the tile 
        x_pos, y_pos - x and y coordinates 
    """
    
    des = record.description
    tile_num = int(des.split(' ')[0].split(':')[4])
    x_pos = int(des.split(' ')[0].split(':')[5])
    y_pos = int(des.split(' ')[0].split(':')[6])
    return tile_num, x_pos, y_pos


def generate_img(x, y, x_min, y_min, x_max, y_max, r, blurred, sigma):

    x_range = x_max - x_min
    y_range = y_max - y_min

    img = np.zeros(shape=(x_range, y_range))
    for i in range(0, len(x)):
        img[min(x_range-1,max(0,int(x[i]-r))):min(x_range-1,max(0,int(x[i]+r+1))),min(y_range-1,max(0,int(y[i]-r))):min(y_range-1,max(0,int(y[i]+r+1)))] = 200

    if blurred:
        img = gaussian(img, sigma=sigma)
    im = Image.fromarray(img)
    new_im = im.convert("L")
#     new_im.save(op_path)
    return new_im

def get_seq_coordinates(fastq_path, tile):
    x_coordinate = []
    y_coordinate = []
    sequences = []
    for record in SeqIO.parse(fastq_path, "fastq"):
        if record is not None:
            tile_num, x_pos, y_pos = get_pos(record)
            seq = str(record.seq)
            if tile_num == tile:
                x_coordinate.append(x_pos/10)
                y_coordinate.append(y_pos/10)
                sequences.append(seq)
#     print('Coordinates are found.')
    return x_coordinate, y_coordinate, sequences
  
def create_FASTQ_image (current_tile, read2, library_seq, seq_threshold):
    """This function is used for creating the image based on FASTQ file. 

    Args:
        current tile(integer): number of the current tile
        read2(boolean):  in case of paired-end sequencing, this flag should be True
        seq(string): template library sequence for selecting library constructs
        seq_threshold(integer): the threshold for a number of matches between a read and the template for selecting library constructs
    Returns:
       Creates the fastq image 
       x_coord - x coordinates from the Fastq file 
       y_coord - y coordinates from the Fastq file
       sequence - list of selected sequences. if read2 =True, it's a 2D array
    """
    fastq_path = fd.askopenfilename(title = "Choose the Read1 FASTQ file")
    
    
    x_coord, y_coord, sequence_1 = get_seq_coordinates(fastq_path,  tile=current_tile)
        
    # Here we rotate the FASTQ coordinated 90deg right and flip them vertically to match the orientation of the FASTQ coordinate
    # system with the smFRET image coordinate system. Please adjust these transformations accordingly
    # (see. Aguire et. al. 2024 DOI: 10.1126/science.adn5371, Fig S16 for the FASTQ axes orientation)
    # x_coord1 = np.subtract(max_y, y_coord)
    # y_coord = np.subtract(max_x, x_coord)
    # x_coord = x_coord1
    
    FQ_coordinates = [x_coord,y_coord]
    if FQ_rotation == 1:
        x_coord = FQ_x_cor*np.array(FQ_coordinates[0])
        y_coord = FQ_y_cor*np.array(FQ_coordinates[1])
        if FQ_x_cor == -1:
            x_coord = np.add(x_coord, max_x)
            
        if FQ_y_cor == -1:
            y_coord = np.add(y_coord, max_y)
            
    if FQ_rotation == -1:
        x_coord = FQ_x_cor*np.array(FQ_coordinates[1])
        y_coord = FQ_y_cor*np.array(FQ_coordinates[0])
        if FQ_x_cor == -1:
            x_coord = np.add(x_coord, max_y)
            
        if FQ_y_cor == -1:
            y_coord = np.add(y_coord, max_x)
            
    
    
    idx = library_index(library_seq, sequence_1, seq_threshold)
    x_coord = [x_coord[i] for i in idx]
    y_coord = [y_coord[i] for i in idx]
    sequence_1 = [sequence_1[i] for i in idx] 
    

    fastq_image = generate_img(y_coord,x_coord, 0, 0, max_x, max_y, 1, True, 1)
    current_direct = fd.askdirectory(title = "Choose the output folder")

    fastq_image.save(os.path.join(current_direct, 'FASTQ_image.png'))
    if read2: 
        fastq_path_2 = fd.askopenfilename(title = "Choose the Read2 FASTQ file") #new
        _, _, sequence_2 = get_seq_coordinates(fastq_path_2,  tile=current_tile)
        sequence_2 = [sequence_2[i] for i in idx]  #new
    
        sequence_1 = [sequence_1, sequence_2]
    return (x_coord, y_coord, sequence_1)
