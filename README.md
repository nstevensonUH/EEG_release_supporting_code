# EEG_release_supporting_code
This repository includes code that supports the public release of the Helsinki 1000 EEG/ECG dataset. This dataset can be found on Zenodo - 

anonymize_dataset_for_github.m - outlines the procedure used to anonymize the dataset, including randomizing filenames, writing generic header files to the EDF files, and quantizing age.
quality_assessment_for_github.m - contains all the code used to calculate the results in the associated publication.

Note this code is developed without sharing in mind. If you are interested in implementing it yourself, please ensure the dataset and associated Matlab code at all in the Matlab path. Use the addpath command in Matlab.

