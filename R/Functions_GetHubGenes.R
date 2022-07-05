## Necessary Libraries
library(hdWGCNA)

## set random seed for reproducibility
getHubGenes_df <- function(csv_file, designation="hub", module="color")  {
  ## Should be dataframe with genes from modules listed with specific set labeled 
  #### make it so function searches for the 'designation' term to find specific column then return respective values
  
  genes <- read.csv(csv_file)
  
  # Search columns for 'designation'
  
  
  # Determine which column contains 'designation' and select rows
  if(any(colnames(genes)==module))  {
    hub_genes <- gene[gene[designation]==designation,]
  }
  else  {
    if (!("hub" %in% unlist(hub_genes))) {
      errorCondition("The designation does not exist in the dataframe, please correct the 'designation' information")
    }
    for (len in 1:ncol(genes))  {
      if(any(genes[,len]==designation)) {
        # Designate correct column with 'designation' values
        hub_genes <- gene[gene[len]==designation,]
      }
    }
  }
  
  # Identify the set of unique 'module' to be analyzed
  if(any(colnames(genes)==module))  {
    module_list <- unique(hub_genes[module])
  }
  else  {
    errorCondition("The module does not exist in the dataframe, please correct the 'module' information")
  }
  
  return(list(hub_genes, module_list))
}


getHubGenes_seurat <- function(seurat_file) {
  ## Dataset should be folder, with Rscripts as a subdirectory containing various projects
  ## Choose one of the projects as the subfolder
  ## WGCNA should have already generated umap_df structure
  
  seurat_obj <- readRDS(file=seurat_file)
  
  umap_df <- GetModuleUMAP(seurat_obj)
  
  if (is.null(umap_df))  {
    errorCondition("WGCNA has not been run on the seurat object or the UMAP does not exist; please run hdWGCNA")
  }
  
  hub_umap_df <- umap_df[umap_df["hub"]=="hub",]
  
  ## Save and return list 
  return(hub_umap_df)
}
