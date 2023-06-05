# -----------------------------------------------------------
# Project: ASD and ADHD - AT project
# This code allows you to brain decode genes i.e. to create a list of genes that are spatially 
# correlated with the imaging maps (e.g. between group t-contrast rfMRI maps) after rigorous
# FDR correction. The statistical modelling simply test the beta values of the gene-imaging 
# correlations of each donor of the Allen Brain Atlas against 0.
# Prior running this analysis, make sure to upload your MRI brain map in NeuroVault 
# (https://neurovault.org/my_collections/?q=). The variable "image" is the NeuroVault ID value.

# Date: June 5, 2023
# -----------------------------------------------------------



# By editing gene_decode.r you can carry out gene decoding either for the whole brain and the cortex.

# edit this, gene_decode.r
source("fMRItogenes_function_gene_decode.r") 

# edit this, your study root folder
dataDir <- getwd()

# edit this, output folder
measure <- "results_of_gene_decoding"

# edit this, the image identifier on NeuroVault (id)
image <- "790618" 

# this is the actual gene decoding calculation
res = gene_decode(dataDir, image, measure, dbs, term, maxterms, prefix, dbs)



