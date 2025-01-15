

def compose (pos):
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
            img_t = np.mean(img_smFRET[green_frames,::], axis = 0)
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
        fig, ax = plt.subplots()
        ax.imshow(combined)
        blobs_log = blob_log(combined, max_sigma=10, num_sigma=10, threshold=300) # Was 1000 for 19/07/2022
    #             Was 300 for 06/09/2022
        CM = []
        r = 3
        [h,w] = red.shape

        for i, blob in enumerate(blobs_log):
            x, y, d = blob
            if x>r and x<(h-r) and y>r and y<(w-r):
                temp = ndimage.measurements.center_of_mass(combined[int(x-r):int(x+r+1),int(y-r):int(y+r+1)])
                CM.append(np.flip(np.add(temp, [x-r,y-r])))

                c = plt.Circle(CM[-1], 3, color="red", linewidth=1, fill=False)
                ax.add_patch(c)
        ax.set_axis_off()
    #     plt.savefig(os.path.join(pos_direct, pos + "_smFRET_peaks.tif"))
        plt.show()
        smFRET_centers = np.array(CM)
    #     if img1.shape[0] == 512: # sometimes slice number is the first dimension, sometimes the third
    #         img1 = np.mean(img1, axis = 2) #averaging by the stack
    #     else:
    #       img1 = np.mean(img1, axis = 0) #averaging by the stack
    # #         img1 = np.mean(img1, axis = 0)
    # #         img1 = np.mean(img1, axis = 2)
    #     img1 = img1.astype("ushort")
    #     img1 = img1[256:,:]
    #     img1 = img1 - si.restoration.rolling_ball(img1, radius=rb_rad)
        smFRET_centers_tr = apriori_tr_original(smFRET_centers)
        smFRET_centers_tr_image_array = np.asarray(generate_img(smFRET_centers_tr[:,1],smFRET_centers_tr[:,0], 0, 0,size1[0,0], size1[0,1], 1, True, 1))
    #     img1 = transform.warp(img1, tr_inv, output_shape = size1[0])
        x1 = [posX[i] for i,x in enumerate(labels) if labels[i]==pos]
        y1 = [posY[i] for i,x in enumerate(labels) if labels[i]==pos]
        id1 = int((maxy-y1[0])/0.34)
        id2 = int((x1[0]-minx)/0.34)

        img[id1:id1+size1[0,0],id2:id2+size1[0,1]] = smFRET_centers_tr_image_array

    img2 = Image.fromarray(img.astype(np.byte))
    img2 = img2.convert('L')
    img2.save(current_direct+'/Combined.png')

    print("DONE")