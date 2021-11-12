# AFM-LineScan-Analyser
MATLAB code to track and analyze objects in AFM Line Scan Kymograph images.

This code was developed in the Scheuring-Lab and used to analyze the motions of GltPH transporters acquired by high-speed AFM line scanning. 

If using this code please refer to and cite:

Matin, T. R., Heath, G. R., Huysmans, G. H., Boudker, O., & Scheuring, S. (2020). Millisecond dynamics of an unlabeled amino acid transporter. Nature communications, 11(1), 1-11.  https://doi.org/10.1038/s41467-020-18811-z 
 
Guide:

Before running the Matlab codes please prepare kymographs as a single image:
Image preparation: Assemble line scanning data into a single tiff image with time in x (left to right), space in y and height as intensity. If the data is saved as a line scan movie this can converted in imageJ using image rotation and montage functions or could be coded in Matlab. 

Using the Matlab codes:
To use these codes, run the following Matlab scripts in editor mode in the following order:
1.	Run S0_LS_Image_analysis to find and track heights of single proteins.

2.	Run S1_1_Kymo_align to remove drift in the image plane (optional)

3.	Run S1_2_track_remover to remove unwanted tracks (optional)

4.	Run S2_LS_Track_Analyzer to analyze number of states from 1 to 3 using a modified StaSI approach.

See notes in each code for specific guidance on how to run each of the above steps. 
