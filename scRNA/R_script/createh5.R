library(scCustomize)

paths = paste0('~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LungM/',
               c('LungM2'))

for (path in paths) {
  
  setwd(path)

  Create_10X_H5(raw_data_file_path = "filtered_feature_bc_matrix", save_file_path = getwd(), save_name = "filtered_feature_bc_matrix")
  Create_10X_H5(raw_data_file_path = "raw_feature_bc_matrix", save_file_path = getwd(), save_name = "raw_feature_bc_matrix")
}