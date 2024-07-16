
# this code carries out the spatial correlation analysis between MRI and genomics data with NeuroVault and outputs:
# - genesample_complete_table.csv --> the complete table with all t-values, p-values, corrected p-values, etc., as in the website 
# - genesample_pos and genesample_neg --> list of genes that are significant at p < 0.05
# - genesample_pos_thres and genesample_neg_thres --> list of genes that are significant at p < 0.05, FDR-corrected

# rather that edit this script, use ‘fMRItogenes_01_decoding_of_brainmaps.r’ to specify your unique inputs


gene_decode <- function (dataDir, image, measure, dbs, term, maxterms, prefix,...){

  # call packages
  require(ggplot2)
  require(dplyr)
  require(rjson)
  require(viridis)
  require(ggpubr)
  
  # do some folder setup
  folderOut <- paste(dataDir,measure,sep="/")
  folderOutTables <- paste(dataDir,measure,"tables",sep="/")
  
  # check if the relevant directories exist, and if not create them
  if (!dir.exists(dataDir))(
    dir.create(dataDir)
  )
  
    if (!dir.exists(folderOut))(
    dir.create(folderOut)
  )
  
  if (!dir.exists(folderOutTables))(
    dir.create(folderOutTables)
  )
  
  # set the working directory
  setwd(dataDir)
  
  # load the gene expression json from Neurovault, using data from the whole brain (/gene_expression/json?mask=full) or limiting the analysis to the cortex (/gene_expression/json?mask=cortex)
  #data = rjson::fromJSON(file=paste("https://neurovault.org/images/",image,"/gene_expression/json?mask=cortex",sep = ""))
  data = rjson::fromJSON(file=paste("https://neurovault.org/images/",image,"/gene_expression/json?mask=full",sep = ""))
  
  # convert output gene list and statistics to dataframe
  df <- data.frame(matrix(t(unlist(data$data)), nrow=length(data$data), byrow=T))
  colnames(df) <- c("symbol","page?","name","t","p","p_corr","var_explained","var_sd")
  
  # make sure the variables have the correct number formatting
  df$t <- as.numeric(as.character(df$t))
  df$p <- as.numeric(as.character(df$p))
  df$p_corr <- as.numeric(as.character(df$p_corr))
  df$var_explained <- as.numeric(as.character(df$var_explained))
  df$var_sd <- as.numeric(as.character(df$var_sd))

  # split into positive and negative gene lists, and threshold at p <0.05 FDR-corrected
  genelist.pos <- df[ which( df$p < 0.05 & df$t >= 0) , ]
  genelist.neg <- df[ which( df$p < 0.05 & df$t <= 0) , ]
  genelist.pos.thres <- df[ which( df$p_corr < 0.05 & df$t >= 0) , ]
  genelist.neg.thres <- df[ which( df$p_corr < 0.05 & df$t <= 0) , ]
  genelist.pos.neg.thres <- df[ which( df$p_corr < 0.05) , ]
  
  # save the resulting tables
  write.table(genelist.pos$symbol,file = paste(folderOut,"/genesample_pos_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg$symbol,file = paste(folderOut,"/genesample_neg_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.pos.thres$symbol,file = paste(folderOut,"/genesample_pos_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg.thres$symbol,file = paste(folderOut,"/genesample_neg_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.pos.neg.thres$symbol,file = paste(folderOut,"/genesample_posORneg_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.csv(df,file = paste(folderOut,"/genesample_complete_table_",measure,".csv",sep = ""), row.names = FALSE)

}
