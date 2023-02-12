library(stringr)
library(Matrix)

### csv data split and merge functions

csvDataSplit = function(data_file, nb_comp, axis = 0) {
  # create data directory
  data_dir = head(strsplit(data_file, "[.]")[[1]], -1)
  dir.create(data_dir, showWarnings = FALSE)
  
  # load data
  df = read.csv(data_file)
  size = ifelse(axis == 0, nrow(df), ncol(df))
  comp_size = ceiling(size/nb_comp)
  
  # subset data
  continue = TRUE; iter = 1
  while (continue) {
    # indexes
    start_index = (iter - 1)*comp_size + 1
    end_index = iter*comp_size
    if (end_index > size) {continue = FALSE; end_index = size}
    idx = start_index:end_index
    
    # save data subset
    if (axis == 0) {
      df_tmp = df[idx, ]
    } else {
      df_tmp = df[, idx]
    }
    file_path = file.path(data_dir, sprintf("comp_%s.csv", iter))
    write.csv(df_tmp, file_path, row.names = FALSE)
    iter = iter + 1
  }
}

# csvDataSplit implementation on large data files for github storage
# 
# url: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248
# data_file = "Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv"
# csvDataSplit(data_file, 22)

csvDataMerge = function(data_dir, axis = 0) {
  df = c()
  print(getwd())
  file_list = list.files(data_dir)
  for (iter in 1:length(file_list)) {
    df_tmp = read.csv( file.path(data_dir, sprintf("comp_%s.csv", iter)) )
    if (axis == 0) {
      df = rbind(df, df_tmp)
    } else {
      df = cbind(df, df_tmp)
    }
  }
  return (df)
}

# csvDataMerge implementation for large data analysis
# 
# data_dir = head(strsplit(data_file, "[.]")[[1]], -1)
# df = csvDataMerge(data_dir)


### dgTMatrix split and merge functions

matrixDataSplit = function(data_file, nb_comp) {
  # create data directory
  data_dir = head(strsplit(data_file, "[.]")[[1]], -1)
  dir.create(data_dir, showWarnings = FALSE)
  
  # load data
  mat = readMM(data_file)
  # df = read.csv(data_file)
  size = length(mat@i)
  comp_size = ceiling(size/nb_comp)
  
  # subset data
  continue = TRUE; iter = 1
  while (continue) {
    # indexes
    start_index = (iter - 1)*comp_size + 1
    end_index = iter*comp_size
    if (end_index > size) {continue = FALSE; end_index = size}
    idx = start_index:end_index
    
    # save data subset
    write.table(mat@i[idx], file = file.path(data_dir, sprintf("i_comp_%s.txt", iter)), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
    write.table(mat@j[idx], file = file.path(data_dir, sprintf("j_comp_%s.txt", iter)), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
    write.table(mat@x[idx], file = file.path(data_dir, sprintf("x_comp_%s.txt", iter)), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
    
    iter = iter + 1
  }
  write.table(mat@Dim, file = file.path(data_dir, "Dim.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
}

# matrixDataSplit implementation on large data files for github storage
# 
# url: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fmatrix%2Emtx%2Egz
# data_file = "GSE113576_matrix.mtx"
# matrixDataSplit(data_file, 10)

matrixDataMerge = function(data_dir) {
  i = c(); j = c(); x = c()
  file_list = list.files(data_dir)
  for (iter in 1:((length(file_list)-1)/3)) {
    i = c(i, scan(file = file.path(data_dir, sprintf("i_comp_%s.txt", iter)), what = integer(), quiet = TRUE))
    j = c(j, scan(file = file.path(data_dir, sprintf("j_comp_%s.txt", iter)), what = integer(), quiet = TRUE))
    x = c(x, scan(file = file.path(data_dir, sprintf("x_comp_%s.txt", iter)), what = double(), quiet = TRUE))
  }
  Dim = scan(file = file.path(data_dir, "Dim.txt"), what = integer(), quiet = TRUE)
  mat = new("dgTMatrix", i = i, j = j, x = x, Dim = Dim)
  return (mat)
}

# matrixDataMerge implementation for large data analysis
# 
# data_dir = head(strsplit(data_file, "[.]")[[1]], -1)
# mat = matrixDataMerge(data_dir)
