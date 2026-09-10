# Cell Movement Analysis

## Project overview

This repository contains data and scripts related to the analysis of cell movement. It is organized in 4 set of different examples that contains the data and scripts for the analysis of different experimental designs, using data obtained with different Fiji analysis (more info in each specific script) for the analysis of movement of _Abeoforma whisleri_ cells. Work done in collaboration with Multicellgenome Lab, CSIC-UPF. 

## Folders and Files

- AtachmentAnalysis:
    - AttachmentAnalysis.R - RScript document containing all code used for analyzing the data obtained with Fiji and several comments explaining the experimental procedure and different pluggins used to obtain and process all data.
    - AttachmentRawData.csv - Data file containing all speed and directionality measurements obtained with Fiji
- BeforeAfterSeparation:
    - BeforeAfterSeparation.R - RScript document containing all code used for analyzing the data obtained with
      Fiji and several comments explaining the experimental procedure and different plugins used to obtain 
      and process all data.
    - Independent Raw Data - Folder containing the .csv data files obtained from Fiji of each Medusa cell analyzed. 
      In here, the Track column have the information about the sample. It is numerated, the order of the numbers
      correspond to: Head Before, Tail Before, Head After, Tail After. The names of the files correspond to the 
      track number of the Medusa cell (example: ResV3.1.2 correspond to the results of Head Before, Tail Before, Head
      After and Tail After of Medusa cell with Track Number 312).
    - Independent Tracks Data - Folder containing all the independent Medusa cell track files. More info in the RScript.
    - MSD Data - Folder containing all the data needed for the MSD analysis. More info in the RScrip. 
    - RawData.csv - Data file containing the concatenation of all files in Independent Raw Data folder for an easy handling.
    - TrackData.csv - Data file containing all the track data from Independent Tracks Data for easy handling.
- CitoskeletonInhibitionAnalysis:
    - CitoskeletonInhibitionAnalysis.R - RScript document containing all code used for analyzing the data obtained with
      Fiji and several comments explaining the experimental procedure and different plugins used to obtain 
      and process all data.
    - TrackData: Folder containing all the data needed for the RScript. Inside this we can find several folders, one for each 
      type of sample used and inside of each one of them we can find all the .csv files containing the data obtained from the
      time-lapses. More information about each one of them in the RScript.
- MedusaMovementAnalysis:
    - CombinedRawData.csv - Data file containing the combination of both, Medusa cells and control data obtained from Fiji. As explained
      in more detail in MedusaMovementAnalysis.R, control and Medusa cells were analyzed independently before comparing both of them. 
      This is the data file needed for the comparison. 
    - CtrRawData.csv - Data file for the data obtained from Fiji for the control cells. 
    - CtrTrackData.csv - Data file containing the track data of the control cells. 
    - MCRawData.csv - Data file for the data obtained from Fiji for the Medusa cells. 
    - MCTrackData.csv -  Data file containing the track data of the Medusa cells. 
    - MedusaMovementAnalysis.R - RScript document containing all code used for analyzing the data obtained with
      Fiji and several comments explaining the experimental procedure and different plugins used to obtain 
      and process all data.

## References

Multicellgenome Lab: https://multicellgenome.com/
