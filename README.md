# MAS
MAS is the Molecularity Analysis Software. This software was developed for high order molecular counting of nanoscopic assemblies with single molecule accuracy from images acquired using confocal, or wide-field, diffraction-limited microscopy. 
## Installation
To install MAS, you will need to install MATLAB 2020a or higher. This version of MATLAB can be downloaded from [here](https://www.mathworks.com/products/matlab.html). Once installed, the following four packages must be installed: signal processing and communications; machine learning and deep learning; maths, statistics and optimization; image processing and computer vision.  
  
Clone this directory in a folder of your choice on your computer.
## Operation
To operate MAS, double-click the `MAS_vx.mlapp` file. The software will load and a GUI will appear.   
  
In the `Data path` panel, press the `Select path` button and a dialog will appear. Browse to the folder containing your data. See below for how to format your data for processing using MAS.  
  
In the `Detection and analysis parameters` panel, flip the `Detect` switch to the `Y` position, insert the number of time slices in the `T (time) slices` field, insert the number of z slices in the `Z (depth) slices` field, insert the number of color channels in the `C (color) slices` field, insert the radius of the Region Of Interest (ROI) containing each assembly (in pixels) in the `ROI radius` field, flip the `Analyze` switch to the `Y` position, insert the stoichiometry of the reference samples (i.e., no. of subunits in each reference assembly) in the `Reference stoichiometry` field, insert the time between each frame in the `Time slice (min)` field and inseet the bin size to be used for Probability Density Function (PDF) underlying the stoichiometry of the unknown assemblies in the `Bin size (mol)` field.
  
Press `Run` button and the text area beneath will show up-to-date information on the status of processing the data. Once the data is processed, the output will be found in the folder named `Analysis` in the folder selected in the `Data path` panel. To repeat, some of the steps above without having to repeat the entire pipeline, flip the switches of the steps not to be repeated to the `N` position.

## Data format
The folder selected in the `Data path` panel should contain two subfolders: a folder named `Calibration` containing z stacks of the reference assemblies and a folder named `Unknown` containing movies (organized as channels, time slices, z slices) of the unknown oligomeric assemblies. All movies should be in an 8-bit multi-tiff format.

## Choice of parameters
Please refer to our original publication below and its supplementary information document for commonly chosen parameters.

## Expected output
A folder named `Analysis` will be exported in the folder selected in the `Data path` panel. This folder will contain 9 files: 
- `Calibration_boxplot.png` file: a picture boxplot of the intensities of the reference assemblies.
- `Calibration_intensity.txt` file: text file of the number of occurances of each intensity value of the reference assemblies.
- `Unknown_molecularity_Xmin.txt` file: text file of the number of occurances of each molecularity value of the unknown assemblies at a certain time point.
- `Unknown_molecularity_PDF_Xmin.txt` file: text files of the probability density function values of the molecularities of the unknown assemblies at a certain time point.
- `Unknown_molecularity_PDF_Xmin.png` file: a picture of the probability density function of the molecularities of the unknown assemblies at a certain time point.
- `Unknown_molecularity_PDF_combined.png` file: a picture of the probability density functions of the molecularities of the unknown assemblies at all time points.
- `Unknown_molecularity_mean/median/emperical_mean/weighted_mean.png` file: a picture of the line plot of the average molecularity (calculated by the mean, median, emperical mean and weighted mean) across time.
- `Unknown_molecularity_mean/median/emperical_mean/weighted_mean.txt` file: text file of the line plot values of the average molecularity (calculated by the mean, median, emperical mean and weighted mean) across time.
- `Unknown_raw.txt` file: text file containing the intensities, positions and times of the unknown assemblies.
## Citing the software
If you use this software in any publication, please cite it as follows:  
