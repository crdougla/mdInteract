####Code by Chris Douglas, PhD UC Irvine ############
###contact crdougla@uci.edu##########################
#####################################################
#####################################################
##Code to evaluate projected drug effectiveness on a module
#####################################################
#####################################################
rm(list=ls()) ##Remove all previous variables / clear the slate


## Load CellMinerCDB Datasets
setwd("C:/Users/Chris/OneDrive/A_UC_LabAdmin_Projects/Research_Projects/R21_PRELIMINARY_GBM_PROJECT_COLLABORATION/Data_Materials/ExternalBioinformaticDatasets/CellMinerCDB/data_NCI-60_exp")
exp <- read.table("data_NCI-60_exp.txt", header=TRUE)
setwd("C:/Users/Chris/OneDrive/A_UC_LabAdmin_Projects/Research_Projects/R21_PRELIMINARY_GBM_PROJECT_COLLABORATION/Data_Materials/ExternalBioinformaticDatasets/CellMinerCDB/drug_NCI-60_act")
drug <- read.table("drug_NCI-60_act.txt", header=TRUE)
 

## Get hub genes
source("C:/Users/Chris/OneDrive/A_UC_LabAdmin_Projects/Research_Projects/R21_PRELIMINARY_GBM_PROJECT_COLLABORATION/Data_Materials/ExternalBioinformaticDatasets/CellMinerCDB/Step01b_GetHubGenes.R")
output <- getHubGenes(Dataset="NatGenetics_2021", subfolder="scWGCNA_4/B", hdWGCNA_file="hdWGCNA_STEP3B.rds")

hub_genes <- output[[1]]

module_list <- output[[2]]

### Determine correlation matrix for gene by drug proliferation inhibition log profile
# Verify which hub genes from modules are in the datasets and initialize module vectors
drug_interactions <- function(hub_genes, module="color", directory)  {
  
  ## wrapped Structure for Input of Module Specific Drug Interactions
  df_module_list <- list()
  
  # Load relevant datasets ######### NEED TO SEE HOW TO REFERENCE WITHIN PACKAGE
  exp <- read.table("data_NCI-60_exp.txt", header=TRUE)
  drug <- read.table("drug_NCI-60_act.txt", header=TRUE)
  
  
  # Loop through all modules
  for (x in unique(hub_genes$module))  {
    
    # Initialize specific interactions for drug with all 'designated' genes
    df_drug <- data.frame(matrix(ncol=10, nrow=nrow(drug)))
    
    # 
    module_genes <- rownames(hub_genes[hub_genes["color"]==x & hub_genes["hub"]=="hub",])
    
    num <- sum(module_genes %in% exp$Gene.name)
    print(paste("Of the 10 hub genes for module: ", x, "; in the drug dataset: ", num, sep=""))
    
    module_exp <- exp[exp$Gene.name %in% module_genes,]
    rownames(module_exp) <- module_exp$Gene.name
    module_exp <- module_exp[,10:ncol(module_exp)]
    
    drug_names <- drug$NAME
    drug_exp <- drug[,5:ncol(drug)]
    
    ## Run correlation of each drug for each gene in module  ########### CAN WE VECTORIZE CORRELATION FUNCTION?
    for (y in 1:length(rownames(module_exp)))  {
      for (z in 1:nrow(drug_exp)) {
        gene <- names(module_exp)[y]
        drug_z <- drug_names[z]
        
        df_drug[z,y] <- cor(unlist(module_exp[y,]), unlist(drug_exp[z,]), use="pairwise.complete.obs")
        
      }
    }
    
    df_drug["average cor"] <- rowMeans(df_drug[,1:nrow(module_exp)])
    df_drug["names"] <- drug_names
    
    df_module_list[[x]] <- df_drug
  }
  
  save(df_module_list, file=paste(directory, "/df_drug_module_list.rda", sep=""))
  
  print(paste("Drug Prediction Data saved in: ", directory, sep=""))
  
  return(df_module_list)
  
}



library(stringi)
library(stringr)
library(sjmisc)
library(ggplot2)
library(cowplot)
library("ggsignif") 

### Take top 40 drugs from each module dataset and apply them
generateDrugBoxPlot <- function (df_module_list, directory)  {
  
  # 
  if (is.null(df_module_list))
  load(paste(directory, "/df_drug_module_list.rda", sep=""))
  dir.create("Drug_activity")
  
}
return_path <- getwd()
load("df_drug_module_list.rda")
dir.create("Drug_activity")
setwd(paste(getwd(), "/Drug_activity", sep=""))
return_path <- getwd()

modules_list <- unique(hub_genes$color)

drug_list <- list()
drug_by_module_id <- list()

for (x in modules_list)  {
  
  ## Create and change to specific module directory
  dir.create(x)
  setwd(paste(getwd(), "/", x, sep=""))
  
  ## Load module specific drug activity dataset and order based on greatest activity - select top 40
  drug_activity_module <- df_module_list[[x]]
  drug_activity_module <- drug_activity_module[order(-drug_activity_module["average cor"]),][1:40,]
  
  module_matrix_by_drug <- data.frame(matrix(ncol=10, nrow=length(modules_list)))
  rownames(module_matrix_by_drug) <- modules_list
  
  # Create list of drugs (e.g. due to empty strings - likely to be less than 40)
  module_x_drugs <- stri_remove_empty(drug_activity_module["names"][1:40,])
  
  drug_vec <- list()
  
  for (y in modules_list)  {
    
    cur_module <- df_module_list[[y]]
    
    # Find the corresponding values for drug activity for top 40 drugs in module x in module y
    drug_vec[[y]] <- cur_module[unlist(cur_module["names"]) %in% module_x_drugs,1:10]

  }
  
  # Nested list, where primary layer refers to module with most n=40 top effective drugs as observed in all modules
  drug_list[[x]] <- drug_vec
  
  # Identities of the top 40 drugs in each module
  drug_by_module_id[[x]] <- module_x_drugs
  
  setwd(return_path)
  
  print(paste("Finished processing drug interactions - module:", x, sep=" "))
}


## Plot Version 3
## Example boxplot
# https://statisticsglobe.com/ggsignif-package-r#:~:text=The%20geom_signif%20function%20also%20enables%20the%20user%20to,line%20size%2C%20and%20text%20size%20as%20shown%20below%3A
# 
# Plot based on each modules top drugs drug_by_module_id[[x]] by referencing each module in 2nd layer of drug_list structure
# Sequential for Error Bars
setwd(return_path)
top <- 4.75
bottom <- 1
err <- seq(from = bottom, to = top, by = (top-bottom)/(length(modules_list)-1))/5 + 0.5

for (x in modules_list)  {
  
  # Top 40 drugs in list
  module_vec <- drug_by_module_id[[x]]
  setwd(x)
  
  # Cycling through top 40 drugs to create list in the respective module subdirectory for module x
  for (z in 1:length(module_vec)) {
    
    drug <- module_vec[z]
    
    # Find module specific response to drug
    response <- vector()
    module_vec_resp <- vector()
    
    for (y in modules_list) {
      value <- as.vector(na.omit(drug_list[[x]][[y]][z,]))
      if (!(is_empty(value)))  {
        response <- append(response, as.numeric(value))
        module_vec_resp <- append(module_vec_resp, unlist(rep(y, length(value))))
      }
    } 
    
    response_df <- as.data.frame(cbind(module_vec_resp, response))
    colnames(response_df) <- c("module", "values")
    
    pairs <- t(combn(unique(response_df$module),2))
    pairs_seq_md <- list()
    for (t in 1:(length(modules_list)-1)) {
      pairs_seq_md[[t]] <- as.character(pairs[t,])
    }
    
    #pairs_md <- split(pairs[pairs[,1]==x,], 1:(length(modules_list)-1))
    
    plot_data<-NULL
    
    png(file=paste(substr(str_replace_all(drug, "[^[:alnum:]]", ""), 1, 30), ".jpeg", sep=""), width = 1000, height = 1000)
    plot_data<-ggplot(data=as.data.frame(response_df), aes(x=factor(module), y=response)) + geom_boxplot() 

    plot_data<- plot_data + geom_signif(comparisons = pairs_seq_md, map_signif_level = TRUE, y_position = err) +
      labs(title = x)
    
    print(plot_data)
    dev.off()
    
  }
  
  print(paste("Generated drug response graphs for module:", x))
  setwd(return_path)
  
}




## TEST sequential list vector for stat comparisons
vec <- seq(1,(length(modules_list)-1), by=1)
pairs_seq <- t(combn(vec,2))
pairs_seq_md <- split(pairs_seq[pairs_seq[,1]==1,], 1:(length(modules_list)-2))

pairs_seq_md <- list()
for (t in 1:(length(modules_list)-1)) {
  pairs_seq_md[[t]] <- as.character(pairs_seq[t,])
}








# Test Sequential for error bar vector
seq(from = 1, to = 4, by = 0.25)




### PLOT VERSION 2
# https://statisticsglobe.com/ggsignif-package-r
library(ggplot2)
library("ggsignif") 

# Sequential for Error Bars
err <- seq(from = 1, to = 4.75, by = 0.25)/5 + 0.5

for (x in drug_bl_gy)  {
  if (str_count(x)<12) {
    print(x)
    jpeg(file=paste(x[1:10], ".jpeg", sep=""), width = 1000, height = 1000)
    ggplot(data=as.data.frame(na.omit(drug_list[[x]])), aes(x=factor(module_vec), y=drug_vec)) + geom_boxplot() + 
      geom_signif(comparisons = list(c("7", "1"), c("7", "2"), c("7", "3"), c("7", "4"), c("7", "5"), c("7", "6"), c("7", "8")
                                     ,c("7", "9"), c("8", "1"), c("8", "2"), c("8", "3"), c("8", "4"), c("8", "5"),
                                     c("8", "6"), c("8", "7"), c("8", "9")), map_signif_level = TRUE, y_position = err)
    dev.off()
    
  }
  
}




### PLOT VERSION 1
for (x in drug_bl_gy)  {
  
  jpeg(file=paste(x[1:10], ".jpeg", sep=""), width = 1000, height = 1000)
  boxplot(drug_vec~module_vec, data=drug_list[[x]])
  dev.off()
  
}


M7_black <- na.exclude(df_mod$black_genes)
num <- sum(M7_black %in% exp$Gene.name)
print(paste("Number of genes in "
))

M7_greenyellow <- na.exclude(df_mod$lightgreen_genes) # name is wrong because I forgot to change back

