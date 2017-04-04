//@File(label="Select the folder with spheroid files", style="directory") sourcedir
//@File(label="Select the folder to save midplane-pictures", style="directory") targetdir
//@File(label="Select the folder to save masks", style="directory") targetdir2

//Copyright (c) 2017 Birga Soetje.
// * 
// * This script is part of the 'spheroid polarity' source code/program.
// *
// * This program is free software: you can redistribute it and/or modify  
// * it under the terms of the GNU General Public License as   
// * published by the Free Software Foundation, version 3.
// *
// * This program is distributed in the hope that it will be useful, but 
// * WITHOUT ANY WARRANTY; without even the implied warranty of 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// * General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program. If not, see <http://www.gnu.org/licenses/>.


//This macro is built to open a folder full of spheroid z-stack images in 4colour fluorescence
//The Order of the fluorescence channels should be 
//   1: gp58(basolateral marker) 
//   2: gp135 (apical marker)
//   3: actin 
//   4: nuclei
//If the channel order differs, an adoption is necesssary!
//It copies the equatorial plane of each spheroid by assuming it in number_of_z/2 +/-3 slices
//It copies the nuclei channel, creates a mask of summed z_slices and determines CenterofMass of the spheroid
//It determines the maximum radius of each spheroid via a projection of all 4 channels, followed by mask creation 
//and calculation of the distance of the outline to the center of mass
//And it determines a potential actin belt via thresholding and creation of a mask, 
//followed by measurement of the position and the area of the actin signal.

//==== global constants =====

//Input: actin_channel - if channelnumber differs please change
var actin_channel = 3;

//Input: nuclei_channel - if channel number differs please change
var nuclei_channel = 4;

// used in the function Actinbelt
// Input a tissue threshold percentage 0-100 
// - here we use 94% to reduce the image to 6% of the signal
var tissueThreshPerc = 94;

//==== user inputs ====
//alternative code for folder selection, replaced by lines 1-3 and 62-64
//Select the source folder - every file in source folder should be an .zvi/.tif or .czi file
//sourcedir = getDirectory("Open the folder with spheroid files");
	
//Choose directory to save midplane pictures and results-file
//targetdir = getDirectory("Select the folder to save midplane-pictures'");	
	
//Choose directory to save nuclei-channel and mask-pictures
//targetdir2 = getDirectory("Select the folder to save masks'");	 

//==== processing ====

sourcedir = sourcedir + File.separator;
targetdir = targetdir + File.separator;
targetdir2 = targetdir2 + File.separator;


spheroids( sourcedir, targetdir, targetdir2 );

// The main function
function spheroids( sourcedir, targetdir, targetdir2 ){
	
	//********* Initialisation *********/
	//Store file list in an array
	filelist = getFileList( sourcedir );
	Array.sort( filelist );
	
	//Prints out all file names within sourcedir
	Array.print( filelist );
	
	//Closes all opened images and clears the results table
	run("Close All");
	run("Clear Results");
	
	//Show file list in results-window
	Array.show("Results", filelist);	

	//Number of spheroid images in the folder
	nSpheroids = lengthOf( filelist );	
	
	//Variable for processed images
	numSpheroids=0;
	
	//Array to store Actin-mask area, coordinates and count of areas (lumen)							
	Actin = newArray( 3 );
	
	//Array to store the CenterofMass, coordinates, area and circularity of the nuclei-projection	
	Center = newArray( 6 );
	
	//Array to store Maximal Radius, CenterofMass, coordinates, 
	//area and circularity of the nuclei-projection															
	MaxRad = newArray( 7 );
	
	//Array to store all results															
	allCenters = newArray();

	//setBatchMode(true);
	//loop over all files and determinatin of all parameters via functions 
	//for Actin, Center and MaximumRadius
	for ( i = 0; i < nSpheroids; i++ ){

		//Check if file is an image file (tif, zvi, czi)
		if ( endsWith( filelist[i], ".zvi") || endsWith( filelist[i], ".tif") || endsWith( filelist[i], ".czi")){		
			
			//Define image name
			currImage = sourcedir + filelist[i];
			print( currImage );
			
			//open image, in case of errors while opening check "Bio-Formats" options and actuality
			run("Bio-Formats (Windowless)", "open=["+ currImage + "]");
			numSpheroids++;
			
			//Arrange Channels to predifined order (see lines 24-27)
			//run("Arrange Channels...", "new=1324");
			
			//Set Unit of Images to pixel / unscaled 
			// - to avoid mixed formats of parameters in case of binned images
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
			ID = getImageID();

			//Function: Determines Actin Parameters: 
			//Area after MaxEntropy Thresholding and Coordinates XM, YM (see below)
			Actin = Actinbelt( ID );		

			//Function: Determines Center / Nuclei Parameters (see below)
			//Parameters XM, YM, Area ,width, height,circ 
			Center = CenterNuclei( ID );

			//Function: Saves the Average Projection of 7 equatorial planes, 
			// determines Center of Mass of whole Spheroid projection: XM spheroid, YM spheroid, 
			// width spheroid, height spheroid, Area spheroid, circ spheroid, maximum radius spheroid
			MaxRad = Midplane_MaxRadius( ID );

			//Concatenating the arrays of actin, nuclei and spheroid paramaters
			Center = Array.concat( Center, MaxRad, Actin);

			//Concatenating the actual results to all results
			allCenters=Array.concat( allCenters,Center );
			Array.show("Backup", filelist, allCenters );
			selectWindow("Backup");
			save( targetdir2 + "ResultsBackup.txt" );
			run("Close");

			//Closes all images
			run("Close All");
		}	
	}	

	//Temporary Output of all Results
	Array.print( allCenters );
	//Create new window for results with title "[Spheroid Parameters]"	
	title1 = "Spheroid_Parameters"; 
	title2 = "[" + title1 + "]";

	//Creates a table with headings containing the results for each spheroid
	//The table content is tab separated and saved under the name "Spheroid_Parameters_Results" as txt document

	run("Table...", "name=" + title2 + " width=350 height=250"); 			
	header = 
		"\\Headings:Image\t" + 
		"X_nuclei\t" + 
		"Y_nuclei\t" + 
		"Area_nuclei\t" + 
		"Width_nuclei\t" + 
		"Height_nuclei\t" + 
		"Circularity_nuclei\t" + 
		"X_spheroid\t" + 
		"Y_spheroid\t" + 
		"Width_spheroid\t" + 
		"Height_spheroid\t" + 
		"Area_spheroid\t" + 
		"Circularity_spheroid\t" + 
		"MaximumRadius_spheroid\t" + 
		"Actin_area\t" + 
		"X_Actin\t" + 
		"Y_Actin\t" + 
		"Count\t" + 
		"Count_after_Watershed\t" + 
		"Circularity_Actin";
	//Parameters of each image
	print( title2, header );  

	//Save each parameter under the name of image within this table
	for (i = 0 ; i < numSpheroids ; i++ ){
		data = filelist[i];
		for ( j = 0; j < 19; j++)
			data = data +"\t" + allCenters[ 19 * i + j ];
		print(title2, data); 
	}
	selectWindow(title1);
	saveAs("Text", targetdir+title1+"_Results.csv");
	 
}


/*
=============== End of the main macro ==============
--------------- Begin of the three functions -----------------
1. CenterNuclei
2. Actinbelt
3. MaxRadius
*/

/**
This function performs a summed z-projection of all slices of the nuclei channel,
creates a black/white mask and determines the parameters (Center of Mass, Area, Dimensions, Shape) 
*/
function CenterNuclei( ID ){
	//Initialisation
	//to choose correct channel for nuclei, get image dimensions and extract image title
	selectImage( ID ); 
	Stack.setChannel( nuclei_channel );
	title = getTitle();
	getDimensions( width, height, channels, slices, frames );
	dotindex = indexOf( title, "." );
	title_prefix = substring( title, 0, dotindex );

	//Duplication of the Nuclei-Channel
	//followed by saving of the image as original-title+"_nuclei.tif"
	run("Duplicate...", "title=[" + title_prefix + "_nuclei.tif] duplicate channels=" + nuclei_channel ); 
	nucleiID = getImageID();							
	save( targetdir2 + title_prefix + "_nuclei.tif" );		

	selectImage( nucleiID );  

	//Z-Projection (Summation) of all slices, without black background
	//and conversion to black/white mask
	run("Z Project...", "projection=[Sum Slices]");
	setOption("BlackBackground", false);
	run("Convert to Mask");						
	maskID = getImageID();
	//Nuclei-Image is closed and mask-image is saved
	if ( isOpen( nucleiID ) ){						
		selectImage( nucleiID );
		close();
	}
	selectImage( maskID );						
	save( targetdir2 + title_prefix + "_mask.tif" );  

	//Set parameters to determine at the created nuclei mask: 
	//area, centroid (Center of Mass), bounding rectangle and shape (value for circularity)
	run("Set Measurements...", "area centroid center bounding shape redirect=None decimal=3"); 

	//If old regions of intrest (ROI) are existing, they are deleted
	if ( roiManager("count") >= 1){  
		roiManager("Deselect");
		roiManager("Deselect");
		roiManager("Delete");
	}
	run("Select None");

	//Determine potential Nuclei-ROIs (size limitation >200 pixel) and add to ROI-Manager
	selectImage( maskID );
	run("Analyze Particles...", "size=200-Infinity pixel add"); 

	//Measure Parameters of the nuclei ROI
	//If more than one nuclei ROI is existing, interpolate parameters over all ROIs
	//Necessary for fragmented or big spheroids
	//Save measured parameters in variables
	run("Clear Results");	
	if ( roiManager("count") > 1 ){  
		count = roiManager("count");
		array = newArray( count );
		for( i = 0; i < count; i++ ) 
			array[i] = i;
		roiManager("Select", array); 
		roiManager("Combine");
		run("Interpolate", "interval=1");
	} else {						
		roiManager("select", 0);
	}
	roiManager("Measure");
	centerX = getResult("XM", nResults-1);
	centerY = getResult("YM", nResults-1);
	Area = getResult("Area", nResults-1);	
	width = getResult("Width", nResults-1);
	height = getResult("Height", nResults-1);	
	circ = getResult("Circ.", nResults-1);		
	
	//Save Parameters in new array "Center" and return this array
	Center = newArray( centerX, centerY, Area, width, height, circ); 

	//Close all images except the original Image
	selectImage( ID ); 
	close("\\Others");
	return Center;
}

/**
This function performs a maximum intensity z-projection of all slices of the actin signal,
creates a black/white mask of it using thresholding at 94% of the signal and 
determines the parameters (Center of Mass, Area) 
*/
function Actinbelt( ID ){
	//Initialisation
	//to choose correct channel for actin, get image dimensions 
	//and extract image title and equatorial plane
	selectImage(ID);
	title = getTitle();
	getDimensions( width, height, channels, slices, frames );
	dotindex =  indexOf( title, "." );
	title_prefix = substring( title, 0, dotindex );

	//Duplicate all slices of the image and project the average intensity
	run("Duplicate...", "title=[" + title_prefix + "_copy.tif] duplicate channels=" + actin_channel );
	run("Z Project...", "projection=[Average Intensity]");
	CopyID = getImageID();

	//Calculate Threshold value for reducing the signal to 6%
	//Save (converted) actinbelt mask in target-folder as original-title+"_actin.tif"
	selectImage( CopyID );

	//8-bit image needed
	run("8-bit"); 
	
	nBins = 256;  
	getHistogram(values, count, nBins);
	size = count.length;
	totalNumberOfPixels = getWidth() * getHeight();

	//Calculate Value, where the cutoff of 94% of the signal is
	tissueValue = totalNumberOfPixels * tissueThreshPerc / 100; 
	print( tissueValue );
	
	// Save cumulative sum of the values within an array
	cumSumValues = count;
	for ( i = 1; i < count.length; i++ )
		cumSumValues[i] += cumSumValues[i-1];
	
	// find position of tissueValue in the array cumSumValues
	for ( i = 1; i < cumSumValues.length; i++ ){
		if ( cumSumValues[i-1] <= tissueValue && tissueValue <= cumSumValues[i]){
			print( i ); // output tissue threshold
			setOption("BlackBackground", true);
			setThreshold( 0, i );
			run("Convert to Mask");
		}
	}
	CopyID = getImageID(); 
	save( targetdir2 + title_prefix + "_actin.tif");		
	selectImage( CopyID );
	run("Invert");
	//delete all old ROIs and analyse for area and center of mass (coordinates)
	if ( roiManager("count") >= 1 ){  
		roiManager("Deselect");
		roiManager("Deselect");
		roiManager("Delete");
	}
	run("Select None");
	run("Set Measurements...", "area centroid center shape median redirect=None decimal=3");
	run("Analyze Particles...", "size=200-Infinity display add");
	count = roiManager("count");

	run("Clear Results");	
	if ( roiManager("count") > 1 ){  
		array=newArray(count);
		for( i = 0; i < count; i++ ) 
			array[i] = i;
		roiManager("Select", array); 
		roiManager("Combine");
		run("Interpolate", "interval=1");
	} else {						
		roiManager("select", 0);
	}
	roiManager("Measure");
	centerX = getResult("XM", nResults-1);
	centerY = getResult("YM", nResults-1);
	Area = getResult("Area", nResults-1);	
	circ = getResult("Circ.", nResults-1);
	
	selectImage( CopyID );
	run("Watershed");
	if ( roiManager("count") >= 1 ){  
		roiManager("Deselect");
		roiManager("Deselect");
		roiManager("Delete");
	}
	run("Select None");
	run("Analyze Particles...", "size=200-Infinity display add");
	count2 = roiManager("count");
	Actin = newArray( Area, centerX, centerY, count, count2, circ );
	
	//Close all images except the original Image
	selectImage( ID );	
	close("\\Others");
	return Actin;

}

//This function performs a projection of the 7 equatorial slices in all 4 channels and saves this image
//Furthermore it projects the sum of all channels and all slices, converts this to a black/white mask 
//and determines the parameters: Area, Centroid(CenterofMass), Bounding Rectangle and Shape
//The Outline of this mask is used to calculate the maximum radius of the spheroid
//As a last step everything outside of the dilated spheroid is set to zero signal via image substraction
function Midplane_MaxRadius( ID ){

	//******* Initialisation ********/

	//extract image title and equatorial plane
	selectImage(ID);
	title = getTitle();
	getDimensions( width, height, channels, slices, frames );
	midplane = round( slices / 2 );
	dotindex = indexOf( title, "." );
	title_prefix = substring( title, 0, dotindex );

	//Duplicate 7 central slices of Image and project the Average Intensity of each channel
	run("Duplicate...", "title=["+title_prefix+"_copy.tif] duplicate channels=1-4 slices="+(midplane-3)+"-"+(midplane+3));
	run("Z Project...", "projection=[Average Intensity]");
	CopyID = getImageID(); 
	
	//for optimal display, loop over all channels to enhance contrast
	for( i = 0; i < channels; i++ ){					
		Stack.setChannel( i + 1 );
		run("Enhance Contrast", "saturated=0.35");
	}
	save( targetdir2 + title_prefix + "_midcopy.tif" );
	CopyID = getImageID();

	//Duplicate all channels of the Image, project everything to one combined plane as Sum
	//and save image as original-title+"_zSum.tif"
	selectImage( ID );
	run("Duplicate...", "duplicate channels=1-4");
	run("Z Project...", "projection=[Max Intensity]");
	run("Grouped Z Project...", "projection=[Sum Slices] group=4");
	ProID = getImageID();
	save( targetdir2 + title_prefix + "_zSum.tif" );

	//Conversion of this projection to a black/white mask
	setAutoThreshold("Mean dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	ID_mask = getImageID();
	save( targetdir2 + title + "_mask.tif");

	//Clearance of the ROI-Manager
	if ( roiManager("count") >= 1 ){  
		roiManager("Deselect");
		roiManager("Deselect");
		roiManager("Delete");
	}
	run("Select None");

	//Setting the parameters to measure: 
	//Area, Centroid(CenterofMass), Bounding Rectangle and Shape
	run("Set Measurements...", "area centroid center bounding shape redirect=None decimal=3");
	
	//Analyse Particles to get spheroid mask
	//followed by inversion and clearance of the outside 
	//to get a black (spheroid) and white (background) mask
	run("Analyze Particles...", "size=500-Infinity pixel display add");

	roiManager("Select", 0);
	setBackgroundColor( 0, 0, 0 ); 
	run("Invert");				
	run("Clear Outside");
	run("Select All");
	run("Invert");

	//Measurement of the spheroid parameters
	roiManager("Select", 0);
	run("Measure");				
	
	//Saving of the measured parameters as variables
	xm1 = getResult("XM", nResults-1); //in pixel
	ym1 = getResult("YM", nResults-1);
	b_w1 = getResult("Width", nResults-1);
	b_h1 = getResult("Height", nResults-1);
	Area = getResult("Area", nResults-1);
	circ = getResult("Circ.", nResults-1);	

	//To determine the maximum radius a new image is generated, which only contains the outline of the spheroid
	//This is realized via flattening of the spheroid-ROI on a white background.
	//After saving of this image the whole image is checked pixel by pixel if this pixel is part of the outline
	//and if so the distance to the center of mass is determined and the largest value for the radius is saved
	radius = 0;

	newImage("Selection", "8-bit white", width, height, 1);
	ID_Sel = getImageID();
	roiManager("Select", 0);
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 1);
	run("Flatten");
	ID_Sel = getImageID();
	selectImage( ID_Sel );
	save( targetdir2 + title_prefix + "_border.tif");
	run("8-bit");

	for( i = 0; i < 360; i++ ){
		//x coordinates
		for( m = 1; m < width; m++ ){ 
			//y coordinates
			for( n = 1; n < height; n++){ 
				if ( getPixel( m, n ) == 0 ){
					//print(m,n);
					rad_temp = 0;
					makeLine( xm1, ym1, m, n );
					run("Measure");
					rad_temp = getResult("Length", nResults-1);
					if ( rad_temp > radius)
						radius = rad_temp;
				}
			}
		}
		radius_a = newArray( xm1, ym1, b_w1, b_h1, Area, circ, radius);

		//The mask image is dilated four times and inverted so that 
		//the slithtly increased spheroid has a value of 0 and the background a maximum of ? in x bit.
		//Afterwards this white/black image is substracted from 
		//the four channels of the average projection of the seven 
		//equatorial slices and by this the background is set to 0.
		selectImage( ID_mask ); 
		if ( roiManager("count") >= 1 ){  
			roiManager("Deselect");
			roiManager("Deselect");
			roiManager("Delete");
		}
		run("Select None");
		// Dilate 4 pixels
		for ( iter = 0; iter < 4; iter++)
			run("Dilate");
		run("Set Measurements...", "area centroid center median redirect=None decimal=3");
		run("Analyze Particles...", "size=500-Infinity pixel display add");

		selectImage( CopyID );
		for( i = 0; i < channels; i++ ){					
			Stack.setChannel( i + 1 );
			roiManager("Select", 0);
			run("Make Inverse");
			run("Measure");
			median = getResult("Median", nResults-1);
			run("Select All");
			run("Subtract...", "value=" + median + " slice");
		}
		MidID = getImageID();
		save( targetdir + title_prefix + ".tif");

		selectImage( ID );
		close("\\Others");

		return radius_a;
	}
}
