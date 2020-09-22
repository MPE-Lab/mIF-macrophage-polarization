#########################################################################################################
#########################################################################################################
#### cell density computation  for (1) M1-like (M1:M2 > threshold), (2) M2-like (M1:M2 < threshold) macrophage cells 
#########################################################################################################
#########################################################################################################
rm(list=ls())

########################################
#### set paths
########################################
root_dir <- file.path(root_path, "test")

input_dir <- file.path(root_dir, 'data')
cell_dir <- file.path(root_dir, 'm1tom2_index')
output_dir <- cell_dir
dir.create(output_dir)

########################################
#### parameters
########################################
megapixelarea <- 0.246# Conversion factor for area (square mm) covered by 1 megapixel at 20x scanning magnification
cell_abundance_cutoff <- 0.3

########################################
#### load cell-level M1:M2 indices (generated from script "M1.M2_cell_density_calculation.r")
########################################
setwd(cell_dir)
M1toM2_indices <- read.csv('M1toM2index_cellLevel.csv', as.is = TRUE)

#### clean up columns
M1toM2_indices <- M1toM2_indices[,c("Sample.Name","Tissue.Category","m1_index","m2_index","m1_to_m2")]

########################################
#### compute cutoff values
########################################
M1_cutoff <- quantile(M1toM2_indices$m1_to_m2, probs=1-cell_abundance_cutoff)
M2_cutoff <- quantile(M1toM2_indices$m1_to_m2, probs=cell_abundance_cutoff)
cutoffs <- c(M2_cutoff,M1_cutoff)

########################################
#### load tissue area data
########################################
setwd(input_dir)
tissue_area <- read.csv('mapping_N_tissueArea.csv', as.is = TRUE)


########################################
#### cell phenotyping based on M1:M2 indices and selected cutoff values
########################################
M1toM2_indices$Phenotype <- ifelse(M1toM2_indices$m1_to_m2 <= cutoffs[1], 'M2', 
                                  ifelse(M1toM2_indices$m1_to_m2 <= cutoffs[2], 'M1M2','M1'))
table(M1toM2_indices$Phenotype)
sum(M1toM2_indices$m1_to_m2 <= cutoffs[1])
sum(M1toM2_indices$m1_to_m2 > cutoffs[2])
nrow(M1toM2_indices)
head(M1toM2_indices)

########################################
#### count cells at core-level
########################################
## unique Sample.Name
sample_nms <- unique(M1toM2_indices$Sample.Name)
 
## loop over inidividual cores/samples
cell_density_perCore <- list()
k <- 1
for (ii in sample_nms){
  #ii <- sample_nms[1]
  cat(ii, '\n')
  data_sub <- M1toM2_indices[M1toM2_indices$Sample.Name == ii, ]
  
  ## ================
  ## count M1-/M2-like cells in tumor (intraepithelial) region
  ## ================
  cellcounttumor <- data.frame('count_M1'= sum(data_sub$Phenotype == 'M1' & data_sub$Tissue.Category=='tumor'),
                               'count_M2'= sum(data_sub$Phenotype == 'M2' & data_sub$Tissue.Category=='tumor'),
                               'count_M1M2'= sum(data_sub$Phenotype == 'M1M2' & data_sub$Tissue.Category=='tumor'))
  colnames(cellcounttumor) <- paste0(colnames(cellcounttumor),'_tumor')
  ## ================
  ## count M1-/M2-like cells in stromal region
  ## ================
  cellcountstroma <- data.frame('count_M1'= sum(data_sub$Phenotype == 'M1' & data_sub$Tissue.Category=='stroma'),
                                'count_M2'= sum(data_sub$Phenotype == 'M2' & data_sub$Tissue.Category=='stroma'),
                                'count_M1M2'= sum(data_sub$Phenotype == 'M1M2' & data_sub$Tissue.Category=='stroma'))
  colnames(cellcountstroma) <- paste0(colnames(cellcountstroma),'_stroma')
  ## ================
  ## count M1-/M2-like cells in overall region
  ## ================
  cellcountall <- data.frame('count_M1'= sum(data_sub$Phenotype == 'M1'),
                             'count_M2'= sum(data_sub$Phenotype == 'M2'),
                             'count_M1M2'= sum(data_sub$Phenotype == 'M1M2'))
  colnames(cellcountall) <- paste0(colnames(cellcountall),'_all')
  
  
  ## ================
  ## get tissue area for current core/sample
  ## ================
  tissue_area_sub <- tissue_area[tissue_area$Sample.Name == ii, ]
  stopifnot(nrow(tissue_area_sub)==1)

  ## ================
  ## compute density = cell counts per tissue region area
  ## ================
  celldensitytumor <- cellcounttumor / (tissue_area_sub$tumorpixels / 1000000 * megapixelarea)
  celldensitystroma <- cellcountstroma / (tissue_area_sub$stromapixels / 1000000 * megapixelarea)
  celldensityall <- cellcountall / ((tissue_area_sub$stromapixels + tissue_area_sub$tumorpixels) / 1000000 * megapixelarea)
  
  ## rename columns
  colnames(celldensitytumor) <- gsub(x=colnames(celldensitytumor), pattern = 'count',replacement = 'density')
  colnames(celldensitystroma) <- gsub(x=colnames(celldensitystroma), pattern = 'count',replacement = 'density')
  colnames(celldensityall) <- gsub(x=colnames(celldensityall), pattern = 'count',replacement = 'density')
  
  ## ================
  ## combine count & density
  ## ================
  count_N_dens <- cbind(cellcounttumor, cellcountstroma, cellcountall,
                        celldensitytumor, celldensitystroma, celldensityall)
  
  cell_density_perCore[[k]] <- count_N_dens
  names(cell_density_perCore)[k] <- tissue_area_sub$core_ids
  k <- k+1
}

ALL_count_N_dens <-  do.call("rbind", cell_density_perCore)
head(ALL_count_N_dens)
setwd(output_dir)
write.csv(ALL_count_N_dens, file = 'coreLvl_count_N_dens.csv', row.names = TRUE)


############################################################
#### concatenate to tumor ids
############################################################
ALL_count_N_dens2 <- merge(tissue_area[,c("core_ids","tumor")], ALL_count_N_dens, by.x = 'core_ids',by.y = 'row.names')
head(ALL_count_N_dens2)
length(unique(ALL_count_N_dens2$tumor))

############################################################
#### calculate tumor-level mean
############################################################
library(dplyr)
ALL_count_N_dens2$core_ids <- NULL
colnames(ALL_count_N_dens2)
tumorAvg <- ALL_count_N_dens2 %>%  group_by(tumor) %>% summarise_all(mean,na.rm=TRUE)
setwd(output_dir)
write.csv(tumorAvg, file = 'tumorLvl_count_N_dens.csv', row.names = FALSE)

