library(tidyverse)

args <- commandArgs(trailingOnly = T)
url_viridian_1 <- as.character(args[1])
url_viridian_2 <- as.character(args[2])
file_path_viridian_1 <- as.character(args[3])
file_path_viridian_2 <- as.character(args[4])
dir_viridian_1 <- as.character(args[5])
dir_viridian_2 <- as.character(args[6])
name_file_descendants <- as.character(args[7])
output_alignment <- as.character(args[8])

ifelse(dir.exists(dir_viridian_1), '', dir.create(dir_viridian_1))
ifelse(dir.exists(dir_viridian_2), '', dir.create(dir_viridian_2))

system(paste0('curl -L -o ', file_path_viridian_1, ' "', url_viridian_1, '"'))
system(paste0('tar -xf ', file_path_viridian_1, ' -C ', dir_viridian_1))
system(paste0('rm ', file_path_viridian_1))
vec_files_viridian_1 <- list.files(paste0(dir_viridian_1, 'Viridian_tree_cons_seqs/'), pattern = "\\.cons\\.fa\\.xz$")
for(file in vec_files_viridian_1){
  system(paste0('unxz -f ', dir_viridian_1, 'Viridian_tree_cons_seqs/', file))
  name_file_unxz <- substr(file, start = 1, stop = nchar(file) - 3)
  system(paste0('seqkit grep -n -f ', name_file_descendants, ' ', dir_viridian_1, 'Viridian_tree_cons_seqs/', name_file_unxz, ' >> ', output_alignment))
  system(paste0('rm ', dir_viridian_1, 'Viridian_tree_cons_seqs/', name_file_unxz))
}

system(paste0('curl -L -o ', file_path_viridian_2, ' "', url_viridian_2, '"'))
system(paste0('tar -xf ', file_path_viridian_2, ' -C ', dir_viridian_2))
system(paste0('rm ', file_path_viridian_2))
vec_files_viridian_2 <- list.files(paste0(dir_viridian_2, 'Viridian_tree_cons_seqs_part_2/'), pattern = "\\.cons\\.fa\\.xz$")
for(file in vec_files_viridian_2){
  system(paste0('unxz -f ', dir_viridian_2, 'Viridian_tree_cons_seqs_part_2/', file))
  name_file_unxz <- substr(file, start = 1, stop = nchar(file) - 3)
  system(paste0('seqkit grep -n -f ', name_file_descendants, ' ', dir_viridian_2, 'Viridian_tree_cons_seqs_part_2/', name_file_unxz, ' >> ', output_alignment))
  system(paste0('rm ', dir_viridian_2, 'Viridian_tree_cons_seqs_part_2/', name_file_unxz))
}
