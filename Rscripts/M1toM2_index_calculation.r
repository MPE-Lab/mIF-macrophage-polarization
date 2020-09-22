rm(list=ls())

########################################
#### parameters
########################################
# 520: CD68 (cytoplasm)
# 540: CD86 (cytoplasm)
# 620: CK (cytoplasm)
# 690: MRC1/CD206 (cytoplasm)
# 570: MAF (nucleus)
# 650: IRF5 (nucleus)
marker_list <- c('Cytoplasm.Opal.520','Cytoplasm.Opal.540','Cytoplasm.Opal.690','Cytoplasm.Opal.620',
                 'Nucleus.Opal.570','Nucleus.Opal.650')

cell_data_filenm <- 'raw_cell_data.csv'

########################################
#### set paths
########################################
root_dir <- getwd()

input_dir <- file.path(root_dir, 'data')
output_dir <- file.path(root_dir,'m1tom2_index')
dir.create(output_dir)

########################################
#### load mapping data: core filenames to tumor_ids
########################################
setwd(input_dir) 
map <- read.csv('mapping_N_tissueArea.csv', row.names = NULL)
map <- map[,c("Sample.Name","tumor")]
head(map)
length(unique(map$tumor))

########################################
#### load raw cell expression data
########################################
setwd(input_dir)
cell_data <- read.csv(file = cell_data_filenm, row.names = NULL)

########################################
#### data cleaning
########################################
## subsetting for cases with patient data
cell_data <- cell_data[cell_data$Sample.Name %in% map$Sample.Name, ]

## subset for macrophage cells 
cell_data <- cell_data[cell_data$cell == 'Macrophage',]

## select cellular compartment 
cell_data<- cell_data[,grep(x=colnames(cell_data) , 
                            pattern = paste0(c("Sample.Name","Tissue.Category",marker_list), collapse = '|'))]
head(cell_data)

## exclude cells with NAs intensity values
cell_data_onlyExprs <- cell_data[,grep(x=colnames(cell_data) , 
                                  pattern = paste0(c(marker_list), collapse = '|'))]
rs <- rowSums(is.na(cell_data_onlyExprs))
sum(rs>0)
cell_data <- cell_data[rs==0,]


########################################
#### on individual markers: expression data scaling to [0 1] range
########################################
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
cell_data_onlyExprs <- cell_data[,grep(x=colnames(cell_data) , 
                                       pattern = paste0(c(marker_list), collapse = '|'))]
scaled_intensity <- apply(cell_data_onlyExprs, MARGIN = 2, range01)

cell_data <- cbind(cell_data[,c("Sample.Name","Tissue.Category")],scaled_intensity)
summary(cell_data)

############################################################
#### compute m1:m2 index = (CD86*IRF5)/(MAF*MRC1)  
############################################################
## compute numerator (CD86*IRF5)
cell_data$m1_index <- cell_data[,"Cytoplasm.Opal.540.Mean..Normalized.Counts..Total.Weighting."] * 
  cell_data[,"Nucleus.Opal.650.Mean..Normalized.Counts..Total.Weighting."]

## compuate denominator (MAF*MRC1)  
cell_data$m2_index <- cell_data[,"Cytoplasm.Opal.690.Mean..Normalized.Counts..Total.Weighting."] * 
  cell_data[,"Nucleus.Opal.570.Mean..Normalized.Counts..Total.Weighting."]

## compuate index
cell_data$m1_to_m2 <- cell_data$m1_index/cell_data$m2_index
## checking for NAs and infinite values
summary(cell_data$m1_to_m2)
sum(is.infinite(cell_data$m1_to_m2 ))

## removing cells with NA or infinite values
cell_data <- cell_data[!is.na(cell_data$m1_to_m2),]
cell_data <- cell_data[!is.infinite(cell_data$m1_to_m2),]

## exclude macrophage in undefined i.e. 'other' tissue region
cell_data <- cell_data[cell_data$Tissue.Category != 'other',]

############################################################
#### append tumor ids
############################################################
cell_data2 <- merge(cell_data, map, by='Sample.Name')
head(cell_data2)
length(unique(cell_data2$tumor))

############################################################
#### saving
############################################################
setwd(output_dir)
write.csv(cell_data2, file = "M1toM2index_cellLevel.csv", row.names = FALSE)

