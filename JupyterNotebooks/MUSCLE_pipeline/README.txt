!NOTE: all the codes were developed for Windows. Running them on other operating systems might require adjustments, e.g. for the different path naming conventions
1. Preparation
In order to run the codes, the following packages should be installed:
os, cv2, tkinter, PIL, skimage, numpy, matplotlib, scipy, cmath, math, sys, json, copy, Bio, sklearn, itertools.
2. Obtaining transformations
Firstly, parameters of a few important transformations should be acquired. 
Transformation which defines the relation between green and red channels. In order to obtain parameters for this transformation the following script can be used: “channel_map_512_2channels_v2.ipynb”. The input for this script is a few movies with fluorescent beads or other fiducial markers that are visible in both Cy3 and Cy5 channels should be used.
!Folder should contain only movies
 After starting the code, the image with 2 channels should appear. Users should click on the same beads from different channels, first click should be at the top channel, second at the bottom channel (e.g. bead1-top, bead1_bottom, bead1_top etc.)
!Don’t click too close to the border 
At least 3 beads need to be manually selected. To finish the clicking step press q on the keyboard. Two binary transformation files will be saved in the same folder as an output (tform_result - forward transformation (from Green to Red), tform_inv - reverse transformation (from Red to Green)). 
A preliminary similarity transformation between smFRET and FASTQ coordinates is required as a starting point for alignment. It can be obtained by calculating the scaling factor based on the FASTQ (0.34um) and smFRET pixel sizes. After the first successful mapping for a given imaging system this transformation can be updated with a similarity transformation generated during the fine alignment stage. Note that the translation part of the transformation should be set to zero, leaving only the rotation and scaling parts.
3. Pipeline usage
Copy the folder “MUSCLE_pipeline” with the library and example (Muscle_data_analysis.ipynb). Note that the latest version of the pipeline can be found at https://github.com/deindllab/MUSCLE.
Here you can find the detailed description of the pipeline, steps’s headers are identical to headers of the cells in Muscle_data_analysis.ipynb: 
! In muscle_sequencing.py find function create_FASTQ_image and change the variable library_seq according to your library sequence.  
0. Import of packages 
 Import all the packages (cell 1 of the example script)

1. Uploading transformations and setting the parameters
Uploading the position list and transformations (described above) can be performed as in the respective cell in the example script.
tform_result - forward transformation (from Green to Red)
tform_inv - reverse transformation (from Red to Green)
Set the important parameters: ALEX (should be true if imaging was performed using alternating laser excitation), read2 (should be true if paired-end sequencing was used), tile (number of the FASTQ tile that will be analyzed)

2. Analysis of the FASTQ file
Generate an image depicting the cluster distribution in a given tile (FASTQ image) and extract x and y coordinates and the sequence of each cluster with the library sequence. 
If sequencing was performed using the paired-end protocol, first choose the Read1 file, then Read2.

!!! check that the flag read2 has the correct value 
x_coord, y_coord, sequence = ma.muscle_sequencing.create_FASTQ_image (current_tile, read2)

more information about the function: help(ma.muscle_sequencing.create_FASTQ_image)

3. Creating the combined image

The coordinates from all FOVs (single-molecule imaging) are used to generate an image depicting the distribution of library molecules across the imaged area.
ma.pks_from_img.combined_image(path_smFRET, POS, ALEX, tr_R2G, apriori_tr_original)
more information about the function: 
help (ma.pks_from_img.combined_image)
ImageJ step
 Combined smFRET image should be aligned to the FASTQ image of a given sequencing tile based on cross-correlation to get an approximate global x and y translation. This step can be accomplished in ImageJ. One way to achieve this is to do the following. First, the two images canvas should be adjusted to the size of the larger of the two images (smFRET composite, Image/Adjust/Canvas size, upper-left). Next, two images can be combined into a single color image (Image/Color/Merge channels). Finally, the two channels can be aligned using a template-matching-ij-plugin (https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin, Align slices in stack). It is easier to align the smFRET composite to the FASTQ image, using a relatively small region of the FASTQ image as a template (100-200um across). The x and y translation that matches the smFRET composite to the FASTQ image need to be provided in the dialog during the next step of the Python pipeline.
4. Alignment
 Here, pairs of clusters and single-molecule coordinates are identified based on spatial proximity and used to generate a refined similarity transformation. The refined transformation is then used to identify final matched pairs. 
rb_rad, r = ma.align.alignment(path_smFRET, POS, ALEX, tr_R2G, apriori_tr, x_coord, y_coord, sequence, read2) #main cycle

At this step the user needs to enter the X and Y displacement values obtained at the ImageJ step.
Also parameters for the background estimation and the half-width of the molecule aperture should be set, we recommend rb_rad = 10, r = 3. These parameters will be used for the trace extraction as well.
5. Trace extraction
SmFRET traces are generated for the matched library molecules using median background estimation:
ma.traces.extract_traces(path_smFRET, read2, rb_rad, r) 
The final output is a set of _traces.mat files containing matched smFRET traces and sequences for each smFRET FOV. Alignment sometimes fails near the edges of the tile. Misaligned FOVs were excluded from further analysis by visual inspection of alignment images in the QC_composites subfolder of the output folder. Further analysis was performed in MatLab 




