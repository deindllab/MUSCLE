// The script assumes the following structure:
// Several folders for different conditions that all contain
// PosNNN subfolders that have red and green subfolders
red_folder = "A_Red2/";
//smFRET_folder = "B_hairpinGreen/"; //Change accordingly!
smFRET_folder_out = "B_Green/"; //Change accordingly!
smFRET_folder = "B_Green_Normal/"; //Change accordingly!
file = "img_000000000.tiff";
ALEX = true; // Flag for choosing the drift correction algorithm

dir = getDirectory("Choose the input directory ");
outdir = getDirectory("Output directory choose you ");
print("dir = ", dir);
list = getFileList(dir);
setBatchMode(true);
//setBatchMode(false);

//print(list.length);

// Removing non-folders from the list in case there are some files alongside folders
for (j=list.length-1; j>=0; j--) {
	if (!endsWith(list[j], "/")) {
		//print("Removed "+list[j]);
		list = Array.deleteIndex(list, j);	
	}
}
Array.sort(list);
positions = getFileList(dir+list[0]);

// Removing non-folders from the list in case there are some files alongside folders
for (j=positions.length-1; j>=0; j--) {
	//print(positions[j]);
	if (!endsWith(positions[j], "/")) {
		//print("removed "+positions[j]);
		positions = Array.deleteIndex(positions, j);
	}
}

//i = 0;
for (i=0; i<positions.length; i++) {
	print("Working on "+positions[i]);
	concat = ""; // Building the string to concatenate all bead images into one
	for (j=0; j<list.length; j++) {
	//for (j=list.length-1; j>=0; j--) {
        filename = dir+list[j]+positions[i]+red_folder+file;
        print(filename);
		open(filename);
		ID = getImageID();
		getDimensions(width, height, channels, slices, frames);
		// Accounting for a possibility of single-frame bead stacks
		if (slices>1) {
			run("Z Project...", "projection=[Average Intensity]");
			rename("beads_"+j); // renaming bead images to reference them during concatenation
			selectImage(ID);
			close();
		} else {
			rename("beads_"+j); // renaming bead images to reference them during concatenation
		}
		

		concat = concat + " image" + (j+1) + "=beads_"+j;
		
		//saveAs(filename1);
		//close();    
	}
	File.makeDirectory(outdir+positions[i]);
	File.makeDirectory(outdir+positions[i]+red_folder);
	File.makeDirectory(outdir+positions[i]+smFRET_folder_out);
	//print(outdir+positions[i]+red_folder);
	run("Concatenate...", concat);
	rename("beads_concat");
	//run("Align slices in stack...", "method=5 windowsizex=430 windowsizey=170 x0=40 y0=300 swindow=0 subpixel=true itpmethod=1 ref.slice=1 show=true");
	run("Align slices in stack...", "method=5 windowsizex=430 windowsizey=180 x0=40 y0=40 swindow=30 subpixel=true itpmethod=1 ref.slice=1 show=true");
	
	//run("Align slices in stack...", "method=5 windowsizex=155 windowsizey=71 x0=190 y0=299 swindow=0 subpixel=true itpmethod=1 ref.slice=1 show=true");
	//run("Align slices in stack...", "method=5 windowsizex=458 windowsizey=130 x0=25 y0=299 swindow=0 subpixel=true itpmethod=1 ref.slice=1 show=true");
	save(outdir+positions[i]+red_folder+file);
	selectImage("beads_concat");
	close();
	numResults = nResults();
	x_shift = newArray(numResults);
	y_shift = newArray(numResults);
	for (j=0; j<numResults; j++) {
		x_shift[j] = getResult("dX",j);
		y_shift[j] = getResult("dY",j);  
	}
	
	concat = " image1 = smFRET_0";
	filename = dir+list[0]+positions[i]+smFRET_folder+file;
	open(filename);
	rename("smFRET_0");
	// Correcting drift within each smFRET movie
	//if (ALEX) {
	//	drift_corr_ALEX();
	//} else {
	//	run("Align slices in stack...", "method=5 windowsizex=430 windowsizey=170 x0=40 y0=300 swindow=5 subpixel=true itpmethod=1 ref.slice=1 show=true");
	//}
	
	for (j=1; j<list.length; j++) {
		filename = dir+list[j]+positions[i]+smFRET_folder+file;
		open(filename);
		rename("smFRET_"+j);
		// Correcting drift within each smFRET movie
		//if (ALEX) {
		//	drift_corr_ALEX();
		//} else {
		//	run("Align slices in stack...", "method=5 windowsizex=430 windowsizey=170 x0=40 y0=300 swindow=5 subpixel=true itpmethod=1 ref.slice=1 show=true");
		//}
		x = x_shift[j-1];
		y = y_shift[j-1];
		concat = concat + " image" + (j+1) + "=smFRET_"+j;
		run("Translate...", "x="+x+" y="+y+" interpolation=Bicubic stack");
	}
	run("Concatenate...", concat);
	//print(concat);
	
	save(outdir+positions[i]+smFRET_folder_out+file);
	close();		
}  
setBatchMode(false);

function drift_corr_ALEX() {
	
	//win = 3; //averaging i-win to i+win
	
	name = getTitle();
	//print(getTitle());
	run("Deinterleave", "how=2 keep");
	//selectWindow(name+" #2");
	selectWindow(name+" #1");
	run("Align slices in stack...", "method=5 windowsizex=430 windowsizey=170 x0=40 y0=300 swindow=10 subpixel=true itpmethod=1 ref.slice=1 show=true");
	selectWindow(name+" #2");
	run("Close");
	selectWindow(name+" #1");
	run("Close");
	len = nResults();
	//print(numResults);
	selectWindow(name);
	x_raw = newArray(len);
	y_raw =newArray(len);
	// Extracting coordinate arrays from lines
	for (k=0; k<len; k++) {
		x_raw[k] = getResult("dX",k);
		y_raw[k] = getResult("dY",k);  
	}
	
	//Preparing coordinate arrays for averaging by extending first and last values
	//temp1 = newArray(win+1);
	//Array.fill(temp1,0);
	//temp2 = newArray(win);
	//Array.fill(temp2,x_raw[len-1]);
	//x_raw = Array.concat(temp1,x_raw,temp2);
	x_raw = Array.concat(0,x_raw);
	//Array.fill(temp1,0);
	//Array.fill(temp2,y_raw[len-1]);
	//y_raw = Array.concat(temp1,y_raw,temp2);
	y_raw = Array.concat(0,y_raw);
	for (k=0; k<len+1; k++) {
//		temp1 = Array.slice(x_raw, i, i + 2*win);	
//		temp2 = Array.slice(y_raw, i, i + 2*win);	
//		Array.getStatistics(temp1, min, max, mean);
		setSlice(2*(k+1));
		run("Translate...", "x="+x_raw[k]+" y="+y_raw[k]+" interpolation=Bicubic slice");    
		print("x="+x_raw[k]+" y=" + y_raw[k]);  
		
	}
	
	for (k=0; k<2*(len+1); k++) {
	    
	}
}
