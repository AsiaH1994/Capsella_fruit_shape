# Investigating changes in capsella fruit shape in response ot increasing temperatures 

### All data for this project can be found in this repository. 

Here we present a site map for this repository: 

The silicula_script_121825.R file contains the R script for statstical analysis and data visulaization.

The files below are used in the R script and the jupyter notebook:

* Temp_fruit_seed_project_010726.csv - information on the 83 inviduals included in this dataset. Data includes seed number, seed weight, later shoot number, by genotype and temperaute.
* genotype_temperature_count.csv - includes the number of individuals planted, the number of individuals that surivied untli senescence, and the ratio of number planted/number survived.
* shape_angles_01072026.csv - necessary information for analysing shape data in the jupyter notebook and measurement data including including PC1/PC2, shape descriptors, and shape measurements
* growth_clean_010726.csv - information on the number of leaves and plant width for each individual over the growing period. Data was collected every 2 days after planting until senescence. 

the flower_test.ipynb jupyter notebook contains code for analyzing Capsella fruit shape

Necessary information for running the shape analysis code includes: 

* image files - included in the image_files folder. All files are zipped. Once unzipped, files should be in jpg format. 
* fruit_all_data.csv - original csv used to analyse all fruit shapes. This csv does not include the subsequent shape descriptor and shape measurement data. This csv is the parent file and was used to create the shape_angles_01072026.csv.

For troubleshooting tips on using the flower_test.ipynb jupyter notebook, please see Hightower, A., S. Hall, R. U. Camacho, A. Papamichail, E. Adamski, C. Colligan, A. Deneen, et al. 2025. Procrustean pseudo-landmark methods in Python to measure massive quantities of leaf shape data. bioRxiv: The Preprint Server for Biology: 2025.08.08.669192. Accepted at APPS October 2025. 



