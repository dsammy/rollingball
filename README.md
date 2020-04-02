# rollingball
Matlab script to correct background of electron backscattered (BSE) images and calculate phase fractions

## Introduction
This repository contains a Matlab script, ```rollingball.m```, that takes as input a set of BSE images and outputs the corresponding thresholded (binary) images. 
Image background is corrected using the rolling ball algorithm introduced by Stanley Sternberg in "Biomedical Image Processing", IEEE Computer, 1983.
Thresholding image segmentation is aided by electron backscattered diffraction (EBSD) phase mapping.


## Citation
The code was written by Mariana C. Mendes Rodrigues. 
Please ```cite``` the following paper if you use it in your own work: M. C. M. Rodrigues, M. Militzer, Application of the rolling ball
algorithm to measure phase volume fraction from backscattered electron images, Materials Characterization 163 (2020) 110273. https://doi.org/10.1016/j.matchar.2020.110273.

For detailed explanations of how to run the code, see below.

## Usage
The following guide assumes that the user has a directory containing the image files. The images should be named using a sequential numerical order. For example, like the following:
```
    - Image1.tif
    - Image2.tif
    - Image3.tif
    - ...
```

Add the BSE image corresponding to the ```EBSD``` scan (with the exact same area) in this same directory and name it as: ```Image_ebsd_map1.tif```, for example. Please note that this BSE image may need to be previously cropped to meet the area of the EBSD map.


### Steps
1. Add the ```rollingball.m``` file to your directory
2. Open Matlab, browse to the root directory and open the ```rollingball.m``` file
3. Edit the root directory of ```myFolder```, in **line 19** and **line 148**
4. Edit the ```tifFilename``` in **line 23** and also in **line 244** from "Ti5553_650C_14min%d.tif" to "Image%d.tif"
5. Change the range value of ```k``` in **line 22** according to your image numbering
6. If necessary, change the ball ```radius``` range in **line 48** (this will depend on the foreground feature size)
7. Edit the output file name in **line 261** from "Binary_image_Ti5553_650C_14min%d.tif" to "Binary_Image%d.tif"
8. Edit the ```tifFilename2``` in **lines 150** and also in **line 208** from "Ti5553_650C_3h_ebsd_map0%d.tif" to "Image_ebsd_map%d.tif"
9. If necessary, change the range value of ```w``` in **line 149** according to your image numbering
10. Execute the script by running ```rollingball``` on the Matlab command line

When the script finishes, your directory will contain the corresponding thresholded images. The calculated phase fractions will appear in a sequential order in the Matlab command line.
