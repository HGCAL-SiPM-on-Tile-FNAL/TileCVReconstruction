# Tile Computer Vision Reconstruction
Code to reconstruct tile edges from multishot and singleshot contour images. The input images are stored in the folder images. They are subfolders to study different cases (e.g. thickness_1 or thickness_2, just different contour thickness). Inside each subfolder, there are three main folders: (1) initialshots, (2) extrashots, and (3) singleshot. 

## Install

Log in to the LPC machines
```
ssh -XY username@cmslpc-sl7.fnal.gov
```
Install python opencv at the LPC 
```
pip install --user opencv-python==4.1.2.30
```

Get code from git
````
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv
git clone https://github.com/HGCAL-SiPM-on-Tile-FNAL/TileReconstruction 
cd TileReconstruction 
````
## Singleshot reconstruction
The script will produce two output folders: plots_singleshot (control plots associated to the filters used) and results_singleshot (fits, shot/contour displaying the fitted corners).
````
python Reconstruction_Singleshot.py
````

## Multishot reconstruction
Make sure to edit the noozle position associated to the shots. The script will produce two output folders: plots_multishot (control plots associated to the filters) and results_multishot (fits, shots/contours displaying the fitted corners).
````
python Reconstruction_Multishot.py
````
