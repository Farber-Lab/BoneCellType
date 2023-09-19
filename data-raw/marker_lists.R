library(usethis)

# cellType_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_cellType.Rdata"))
# encode greek letters to UTF-8
cellType_list = lapply(cellType_list, enc2utf8)
usethis::use_data(cellType_list, overwrite = TRUE)

# cellTypeSource_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_cellType_and_sourceDB.Rdata"))
# encode greek letters to UTF-8
cellTypeSource_list = lapply(cellTypeSource_list, enc2utf8)
usethis::use_data(cellTypeSource_list, overwrite = TRUE)

# tissueSource_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_tissue_and_sourceDB.Rdata"))
# encode greek letters to UTF-8
tissueSource_list = lapply(tissueSource_list, enc2utf8)
usethis::use_data(tissueSource_list, overwrite = TRUE)

# tissue_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_tissue.Rdata"))
# encode greek letters to UTF-8
tissue_list = lapply(tissue_list, enc2utf8)
usethis::use_data(tissue_list, overwrite = TRUE)
