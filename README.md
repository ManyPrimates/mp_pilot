# ManyPrimates Pilot

Project repository for the ManyPrimates Pilot study. For more information about ManyPriamtes visit: https://sites.google.com/view/manyprimates/home

## Note

Datafiles and code are constantly being updated as new data files are submitted. This message will be deleted as soon as all data files have been submitted.

## Structure

Raw data files, as submitted from each lab, are stored in data/raw_data/. The merged dataset is stored in data/merged_data. The code for data pre-processing prior to merging the datasets can be found in analysis/01_data_processing.Rmd. In this files, we also insert life expectancy for each species and calculate the normed age per subject and normalize this variable within each species. The file analysis/02_visualization.Rmd contains code for some preliminary visualizations.
