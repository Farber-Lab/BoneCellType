library(usethis)
markerDatabase = read.delim("https://raw.githubusercontent.com/kkupkova/Mouse-bone-markers/main/BONE_DATABASE.tsv")
markerDatabase = markerDatabase[, !(names(markerDatabase) == "type")]
usethis::use_data(markerDatabase, overwrite = TRUE)
