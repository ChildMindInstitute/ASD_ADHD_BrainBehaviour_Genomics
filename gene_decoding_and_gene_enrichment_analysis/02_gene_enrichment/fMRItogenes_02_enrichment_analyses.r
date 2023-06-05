# -----------------------------------------------------------
# Project: ASD and ADHD - AT project
# this code tests significant overlapping between the list previously obtained with brain decoding analysis and
# the list resulting from genes that are found to be differentialy expressed in the
# autism cortex.
# Statistical testing of overlapping between the two list genes is performed with hypergeometric testing.  

# Date: June 5, 2023
# -----------------------------------------------------------

# background totals to use for enrichment analyses 
ns_gex_background = 20787

# source main enrichment script, edit /path/genelistOverlap.r
source("fMRItogenes_function_genelistOverlap.r")

# edit this, this is the path of the root folder
rootpath = "gene_decoding_and_gene_enrichment_analysis"

# list to modify, e.g. Satterstrom ASD+ADHD genes combining PTV+Mis (Cases no ID vs Controls, via Fisher's exact test)
asd_adhd_mut = unique(as.character(read.delim(file.path(rootpath,"Satterstrom_genelist_alterations_ASD_ADHD.txt"),header = FALSE)$V1))

# this reads the brain decoded gene list, if needed edit the folder name (path) and the list name (txt file) 
decoded_gene_list_positive_fdr05 = unique(as.character(read.delim(file.path(rootpath, "results_of_gene_decoding", "genesample_pos_thres_results_of_gene_decoding.txt"), header = FALSE)$V1))
decoded_gene_list_negative_fdr05 = unique(as.character(read.delim(file.path(rootpath, "results_of_gene_decoding", "genesample_neg_thres_results_of_gene_decoding.txt"), header = FALSE)$V1))
decoded_gene_list_positivenegative_fdr05 = unique(as.character(read.delim(file.path(rootpath, "results_of_gene_decoding", "genesample_posORneg_thres_results_of_gene_decoding.txt"), header = FALSE)$V1))

# this output odds ratio, p-value, number of overlapping genes and other results of the enrichment analysis 
# between ASD DE co-expression modules and the brain decoded gene list
results_enrichment_positive = genelistOverlap(decoded_gene_list_positive_fdr05, asd_adhd_mut, ns_gex_background)
results_enrichment_negative = genelistOverlap(decoded_gene_list_negative_fdr05, asd_adhd_mut, ns_gex_background)
results_enrichment_positivenegative = genelistOverlap(decoded_gene_list_positivenegative_fdr05, asd_adhd_mut, ns_gex_background)


# this writes odds ratio, p-value, number of overlapping genes and other results of the enrichment analysis 
# between ASD DE co-expression modules and the brain decoded gene list
sink("results_enrichment_positive.txt")
print(results_enrichment_positive)
sink()


sink("results_enrichment_negative.txt")
print(results_enrichment_negative)
sink()

sink("results_enrichment_positivenegative.txt")
print(results_enrichment_positivenegative)
sink()


# this outputs the list of enriched genes between ASD DE co-expression modules and the decoded gene list
genelist_enrichment_positive = sort(results_enrichment_positive[[1]]$overlapping_genes)
genelist_enrichment_negative = sort(results_enrichment_negative[[1]]$overlapping_genes)
genelist_enrichment_positivenegative = sort(results_enrichment_positivenegative[[1]]$overlapping_genes)


# this writes the list of enriched genes between ASD and/or ADHD altered genes and the decoded gene list to txt file
write(genelist_enrichment_positive, file = "genelist_enriched_positive.txt")
write(genelist_enrichment_negative, file = "genelist_enriched_negative.txt")
write(genelist_enrichment_positivenegative, file = "genelist_enriched_positivenegative.txt")


