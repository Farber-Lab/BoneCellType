library(usethis)

# cellType_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_cellType.Rdata"))
usethis::use_data(cellType_list, overwrite = TRUE)

# cellTypeSource_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_cellType_and_sourceDB.Rdata"))
usethis::use_data(cellTypeSource_list, overwrite = TRUE)

# tissueSource_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_tissue_and_sourceDB.Rdata"))
usethis::use_data(tissueSource_list, overwrite = TRUE)

# tissue_list
load(url("https://github.com/kkupkova/Mouse-bone-markers/raw/main/BONE_DATABASE_lists/bone_tissue.Rdata"))
usethis::use_data(tissue_list, overwrite = TRUE)
