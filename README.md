# ManyPrimates Pilot

Project repository for the ManyPrimates Pilot study. For more information about ManyPriamtes visit: https://sites.google.com/view/manyprimates/home


## Structure

Raw data files, as submitted from each lab, are stored in data/raw_data/. 

The merged dataset is stored in data/merged_data. The code for data pre-processing prior to merging the datasets can be found in analysis/01_data_processing.Rmd. In this files, we also insert life expectancy for each species and calculate the normed age per subject. 

The file analysis/02_visualization.Rmd contains code for some preliminary visualizations. 

In analysis/04_inference.Rmd we are running simple inferential statistics. 

The code for the model we use to analyse the combined data set is in analysis/05_mixedmodel.Rmd. 
