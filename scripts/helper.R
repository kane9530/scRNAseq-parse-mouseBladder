# color palette taken from archR package. Shifted to start from 0.
stallion = c("0"="#D51F26","1"="#272E6A","2"="#208A42","3"="#89288F","4"="#F47D2B", "5"="#FEE500","6"="#8A9FD1","7"="#C06CAB","18"="#E6C2DC",
             "8"="#90D5E4", "10"="#89C75F","11"="#F37B7D","12"="#9983BD","13"="#D24B27","14"="#3BBCA8", "15"="#6E4B9E","16"="#0C727C", "17"="#7E1416","9"="#D8A767","19"="#3D3D3D")

SaveFigure <- function(res_path, plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(res_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(res_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(rds_path, object, name){
  saveRDS(object, paste0(rds_path, name, ".RDS"))
}

ReadObject <- function(rds_path,name){
  readRDS(paste0(rds_path, name, ".RDS"))
}

# Function to change x labels of individual violin plots
VlnPlot_2 <- function(object, features.plot, title) {
  
  # Main function
  main_function <- function(object = object, features.plot = features.plot) {
    VlnPlot(object = object, features = features.plot) + 
      theme(legend.position="none") +
      theme(axis.text.x=element_blank())
  }
  
  # Apply main function on all features
  p <- lapply(X = features.plot, object = object,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(features.plot))))
}

loadIntoSeurat <- function(path_to_sample){
    
  # Input: Path to the DGE_filtered folder corresponding to the 
  # sample of interest. The folder should contain a DGE.mtx file, a cell_metadata.csv 
  # file and an all_genes.csv file.
  # 
  # Output: A seurat object with the appropriate metadata added to seuratObj@meta.data
  # entry
    
  res_mat <- ReadParseBio(paste0(path_to_sample))
  print(paste("[loadIntoSeurat] -  Loaded parseBio count matrix from", path_to_sample))
  res_mat <- as(res_mat, "dgCMatrix")
  genes <- genes <- read.csv(paste0(path_to_sample, "all_genes.csv"), header=TRUE) %>% dplyr::pull("gene_name")
  print("[loadIntoSeurat] - Loaded parseBio features!")
  cell_metadata <- read.csv(paste0(path_to_sample, "cell_metadata.csv"))
  cell_metadata <- cell_metadata %>% tidyr::separate(sample, c("Condition", "Type"))
  print("[loadIntoSeurat] - Loaded parseBio cell metadata!")
  
  # Setting genes as rownames for count matrix
  stopifnot(dim(res_mat)[1] == length(genes))
  rownames(res_mat) <- genes
  table(rownames(res_mat) == "")
  # Check to see if empty gene names are present, add name if so.
  rownames(res_mat)[rownames(res_mat) == ""] <- "unknown"
  print("[loadIntoSeurat] - Assigned genes as rownames in count matrix!")
  
  # Including sublibrary identity into metadata
  sublibrary_ident <- sub(".*/(Sublibrary[1|2])/.*","\\1",path_to_sample)
  cell_metadata$Sublibrary <- sublibrary_ident
  print(sublibrary_ident)
  print("[loadIntoSeurat] - Included sublibrary identity into metadata!")
  
  # Set cell metadata rownames to colnames of count matrix
  rownames(cell_metadata) <- cell_metadata$bc_wells
  print("[loadIntoSeurat] - Set metadata rownames to cell names in count matrix!")
  # Identify Sample identity
  sample_ident <- sub(".*/(KO.*|WT.*)/DGE.*","\\1",path_to_sample)
  # Create Seurat object 
  seuratObj <- CreateSeuratObject(counts = res_mat, 
                                  meta.data = cell_metadata,
                                  row.names= rownames(res_mat))
  # Set identity of each cell to sublib+sampleidentity
  Idents(seuratObj) <- paste(sublibrary_ident,sample_ident, sep="_")
  print("[loadIntoSeurat] - Created seurat object!")
  return(seuratObj)
}

generateQCmetrics <- function(seuratObj, res_path){
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt")
  seuratObj[["percent.ribosomal"]] <- PercentageFeatureSet(seuratObj, pattern = "^Rp[sl]")
  
  qc.metrics <- as_tibble(
    seuratObj[[]],
    rownames="Cell.Barcode"
  ) 
  
  ## nCount_RNA
  ncountRNAhist <- qc.metrics %>%
    ggplot(aes(log10(nCount_RNA))) + 
    geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
    ggtitle("Distribution of UMI counts") +
    theme_minimal()+
    theme(plot.title = element_text(size = 12))
  
  ## nfeatures_RNA
  nfeatureshist <- qc.metrics %>%
    ggplot(aes(log10(nFeature_RNA))) + 
    geom_histogram(binwidth = 0.1, fill="yellow", colour="black") +
    ggtitle("Distribution of number of detected genes") +
    theme_minimal()+
    theme(plot.title = element_text(size = 12))
  nfeatureshist
  
  ## percent.ribo
  ncountRibohist <- qc.metrics %>%
    ggplot(aes(percent.ribosomal)) + 
    geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage ribosomal genes") +
    theme_minimal()+
    theme(plot.title = element_text(size = 12))
  
  ## percent.mt
  percent.mt_hist <- qc.metrics %>%
    ggplot(aes(percent.mt)) + 
    geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage Mitochondrial genes") +
    theme_minimal()+
    theme(plot.title = element_text(size = 12))
  percent.mt_hist
  
  combined_qc <- cowplot::plot_grid(ncountRNAhist, nfeatureshist, ncountRibohist, percent.mt_hist) 
  #print(unique(Idents(seuratObj)))
  #png(paste0(res_path, unique(Idents(seuratObj)), "_QC_plots",".png"),
  #    width = 12, height = 6, units = "in", res = 200)
  #combined_qc
  SaveFigure(res_path, combined_qc, paste(unique(Idents(seuratObj)),"QC_plots",sep="_"), width = 12, height = 6)
  print("[generateQCmetrics] - Output combined figure!")
  return(seuratObj)
}
