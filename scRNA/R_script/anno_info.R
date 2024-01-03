

path = '~/Project/MultiOmics/data/snRNA/Object/'
sample = 'summary'

create_matrix = function(path, sample) {
  predictions_BCatlas_major <-
    read.csv(file.path(path, sample, "annotation/predictions_BCatlas_major.csv"),
             row.names = 1)
  predictions_BCatlas_minor <-
    read.csv(file.path(path, sample, "annotation/predictions_BCatlas_minor.csv"),
             row.names = 1)
  predictions_BCatlas_subset <-
    read.csv(
      file.path(path, sample, "annotation/predictions_BCatlas_subset.csv"),
      row.names = 1
    )
  
  predictions_HBCA_author_cell_type <-
    read.csv(
      file.path(
        path,
        sample,
        "annotation/predictions_HBCA_author_cell_type.csv"
      ),
      row.names = 1
    )
  predictions_HBCA_cell_type <-
    read.csv(
      file.path(path, sample, "annotation/predictions_HBCA_cell_type.csv"),
      row.names = 1
    )
  predictions_HBCA_broad_cell_type <-
    read.csv(
      file.path(
        path,
        sample,
        "annotation/predictions_HBCA_broad_cell_type.csv"
      ),
      row.names = 1
    )
  
  scanvi_BCatlas <-
    read.csv(file.path(path, sample, "annotation/scanvi_BCatlas.csv"),
             row.names = 1)
  manual_annotation <-
    read.csv(file.path(path, sample, "annotation/manual_annotation.csv"),
             row.names = 1)
  
  anno.matrix = cbind(
    BCatlas_major = predictions_BCatlas_major,
    BCatlas_minor = predictions_BCatlas_minor,
    BCatlas_subset = predictions_BCatlas_subset,
    HBCA_broad_cell_type = predictions_HBCA_broad_cell_type,
    HBCA_cell_type = predictions_HBCA_cell_type,
    HBCA_author_cell_type = predictions_HBCA_author_cell_type,
    scanvi_BCatlas = scanvi_BCatlas,
    manual_annotation = manual_annotation
  )
  
  anno.matrix = anno.matrix[colnames(anno.matrix) %>% grep('clustering', ., invert = T, value = T) %>% grep('predicted', ., invert = T, value = T)]
  return(anno.matrix)
  
}

anno.matrix = create_matrix(path, sample)

library(gtsummary)
anno.matrix %>% tbl_summary()

#table(anno.matrix$manual_celltype_annotation,anno.matrix$HBCA_broad_cell_type.majority_voting)
