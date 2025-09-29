shared_allen_anno <-
  read.csv(file = "data-raw/AllenCCFv3_anno.txt", header = T, sep = ":")

usethis::use_data(shared_allen_anno, overwrite = T)
