// The script assumes the following structure:
// A folder that contains PosNNN subfolders that have green subfolders
// The goal is to extract red excitation frames and save them as a red subfolder
red_folder = "A_Red2/";
//smFRET_folder = "B_hairpinGreen/"; //Change accordingly!
smFRET_folder = "B_Green/"; //Change accordingly!
file = "img_000000000.tiff";

dir = getDirectory("Choose the input directory ");
//outdir = getDirectory("Output directory choose you ");
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
positions = list;

//i = 0;
for (i=0; i<positions.length; i++) {
	print("Working on "+positions[i]);
//	concat = ""; // Building the string to concatenate all bead images into one
//	for (j=0; j<list.length; j++) {
	//for (j=list.length-1; j>=0; j--) {
    filename = dir+positions[i]+smFRET_folder+file;
    //print(filename);
	open(filename);
	run("Deinterleave", "how=2");
	close();

	//ID = getImageID();
	
	File.makeDirectory(dir+positions[i]+red_folder);
	save(dir+positions[i]+red_folder+file);
	close();
	//getDimensions(width, height, channels, slices, frames);
	// Accounting for a possibility of single-frame bead stacks
}

setBatchMode(false);
