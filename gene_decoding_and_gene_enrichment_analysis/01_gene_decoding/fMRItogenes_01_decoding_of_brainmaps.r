# -----------------------------------------------------------
# Project: ASD and ADHD - AT project
# This code allows you to do gene decoding of brainmaps, i.e. to create a list of genes that are spatially 
# correlated with the neuroimaging maps (e.g., group t-contrast resting state functional MRI maps) after rigorous
# FDR correction. The statistical modelling simply tests the beta values of the gene-imaging 
# correlations of each donor of the Allen Brain Atlas against 0.
# Prior running this analysis, make sure to upload your MRI brain map in NeuroVault 
# (https://neurovault.org/my_collections/?q=). The variable "image" is the NeuroVault ID value.

# Date: June 5, 2023
# -----------------------------------------------------------



# By editing fMRItogenes_function_gene_decode.r you can carry out gene decoding either for the whole brain or the cortex. 
# Currently, it is set to carry out decoding for the whole brain. 

# edit this to give the path to the fMRItogenes_function_gene_decode.r file
source("fMRItogenes_function_gene_decode.r") 

# edit this to your study root folder
dataDir <- getwd()

# edit this to the name of your desired output folder
measure <- "results_of_gene_decoding"

# edit this to the image identifier on NeuroVault (id)
image <- "790618" 

# this line runs the actual gene decoding calculation
res = gene_decode(dataDir, image, measure, dbs, term, maxterms, prefix, dbs)



