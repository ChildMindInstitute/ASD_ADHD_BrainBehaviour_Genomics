
# this code carries out the spatial correlation analysis between MRI and genomics data with NeuroVault and outputs:
# - genesample_complete_table.csv --> the complete table with all t- and p-values and the rest, as in the website 
# - genesample_pos and genesample_neg --> list of genes that are significant at p < 0.05
# - genesample_pos_thres and genesample_neg_thres --> list of genes that are significant at p < 0.05, FDR-corrected

# you don't need to edit this code


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
  
  # check if it exists and if not make it so
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
  
  # load the json from Neurosynth, use data from the whole brain or limit the analysis to the cortex
  #data = rjson::fromJSON(file=paste("https://neurovault.org/images/",image,"/gene_expression/json?mask=cortex",sep = ""))
  data = rjson::fromJSON(file=paste("https://neurovault.org/images/",image,"/gene_expression/json?mask=full",sep = ""))
  
  # I don't like lists so convert to a usable dataframe (there's probably a better way to do this...)
  df <- data.frame(matrix(t(unlist(data$data)), nrow=length(data$data), byrow=T))
  colnames(df) <- c("symbol","page?","name","t","p","p_corr","var_explained","var_sd")
  
  # now make sure they have the correct format again
  df$t <- as.numeric(as.character(df$t))
  df$p <- as.numeric(as.character(df$p))
  df$p_corr <- as.numeric(as.character(df$p_corr))
  df$var_explained <- as.numeric(as.character(df$var_explained))
  df$var_sd <- as.numeric(as.character(df$var_sd))

  # split positive and negative and threshold
  genelist.pos <- df[ which( df$p < 0.05 & df$t >= 0) , ]
  genelist.neg <- df[ which( df$p < 0.05 & df$t <= 0) , ]
  genelist.pos.thres <- df[ which( df$p_corr < 0.05 & df$t >= 0) , ]
  genelist.neg.thres <- df[ which( df$p_corr < 0.05 & df$t <= 0) , ]
  genelist.pos.neg.thres <- df[ which( df$p_corr < 0.05) , ]
  
  # save the tables
  write.table(genelist.pos$symbol,file = paste(folderOut,"/genesample_pos_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg$symbol,file = paste(folderOut,"/genesample_neg_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.pos.thres$symbol,file = paste(folderOut,"/genesample_pos_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg.thres$symbol,file = paste(folderOut,"/genesample_neg_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.pos.neg.thres$symbol,file = paste(folderOut,"/genesample_posORneg_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.csv(df,file = paste(folderOut,"/genesample_complete_table_",measure,".csv",sep = ""), row.names = FALSE)

}
