# SAS
SAS is the Stoichiometry Analysis Software. This software was developed for high order subunit counting of biomolecular complexes at the single molecule level acquired using wide-field diffraction-limited microscopy. 
## Installation
To install SAS, you will need to install MATLAB 2020a or higher. This version of MATLAB can be downloaded from [here](https://www.mathworks.com/products/matlab.html). Once installed, the following four packages must be installed: signal processing and communications; machine learning and deep learning; maths, statistics and optimization; image processing and computer vision.  
  
Clone this directory in a folder of your choice on your computer.
## Operation
To operate SAS, double-click the `SAS_vx.mlapp` file. The software will load and a GUI will appear.   
  
In the `Data path` panel, press the `Select path` button and a dialog will appear. Browse to the folder containing your data. See below for how to format your data for processing using SAS.  
  
In the `Detection parameters` panel, flip the `Detect` switch to the `Y` position, insert the pixel size of the detection camera (in nanometers) in the `Camera pixel size` field, insert the quantum efficiency of the detection camera (as percentage) at the emission wavelength of the used fluorophore in the `Camera quantum efficiency` field, insert the offset of the detection camera (in signal units) in the `Camera offset` field, insert the electron multiplication gain of the detection camera (in in electrons / photons) in the `Camera EM gain` field, and a maximum value of the standard deviation in the point spreaf function of the detected particles (in nanometers) in the `Maximum sigma` field. All camera parameters can be obtained from the supplier. 
  
In the `Processing parameters` panel, flip the `Process` switch to the `Y` position and insert the radius of the Region Of Interest (ROI) containing each particle (in pixels) in the `ROI radius` field. 
  
In the `Annotation parameters` panel, flip the `Annotate` switch to the `Y` position and insert the minimum height of the normalized gradient traces (as percentage) in the `Minimum peak height` field.
  
In the `Analysis parameters` panel, flip the `Analyze` switch to the `Y` position, insert the efficiency of fluorophore labeling (as percentage) in the `Labeling efficiency` field and the maximum number of gaussians used for fitting the unknown species in the `Maximum Gaussian mixtures` field.  
  
Press `Run` button and the text area beneath will show up-to-date information on the status of processing the data. Once the data is processed, the output will be found in the folder named `Analysis` in the folder selected in the `Data path` panel. To repeat, some of the steps above without having to repeat the entire pipeline, flip the switches of the steps not to be repeart to the `N` position.

## Data format
The folder selected in the `Data path` panel should contain two subfolders: a folder named `Calibration` containing all calibration movies (ideally containing monomeric particles) and a folder named `Unknown` containing movies of the unknown oligomeric species. All movies should be in an 8-bit multi-tiff format.

## Choice of parameters
Please refer to our original publication below for how to best choose your parameters. The `Minimum peak height` parameters has to be adjusted once by processing the data and checking the obtained calibration data as described below.

## Expected output
A folder named `Analysis` will be exported in the folder selected in the `Data path` panel. This folder will contain 9 files: 
- `Calibration_PDF.png` file: a picture plot of the fitted calibration Probability Density Function (PDF). The curve should ideally show a single peak corresponding to the monomeric distribution. If more than a single peak is shown, or if the fitted curve does not fit the distribution appropriately, adjust the `Minimum peak height` field as described above.
- `Calibration_PDF.txt` file: text file containing the number data of the `Calibration_PDF.png` picture plot for reproduction.
- `Calibration_raw.txt` file: text file containing the intensities and positions of the particles selected to construct the calibration PDF.
- `Unknown_PDF.png` file: a pixture plot of the fitted Probability Density Function (PDF) of the unknown oligomeric species.
- `Unknown_PDF.txt` file: text file containing the number data of the `Unknown_PDF.png` picture plot for reproduction.
- `Unknown_raw.txt` file: text file containing the intensities and positions of the particles used in constructing the PDF of the unknown oligomeric species.
- `Unknown_species proportion(before labelling correction).png` file: a pixture plot of the bar chart of the stoichiometry of the processed complexes before correction of labeling efficiency.
- `Unknown_species proportion(after labelling correction).png` file: a pixture plot of the bar chart of the stoichiometry of the processed complexes after correction of labeling efficiency.
- `Unknown_species proportion(summary).txt` file: text file containing the number data of the `Unknown_species proportion(before labelling correction).png` and `Unknown_species proportion(after labelling correction).png` picture plots.
## Citing the software
If you use this software in any publication, please cite it as follows:  
  
Danial, JSH, Quintana, Y, Ros, U, Shalaby, R, Margheritis, EG, Chumpen, S, Ungermann, C, Garcia-Saez, AJ, Cosentino, K. Systematic assessment of the accuracy of subunit counting in biomolecular complexes using automated single molecule brightness analysis.
