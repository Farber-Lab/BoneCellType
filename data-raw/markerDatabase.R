library(usethis)
markerDatabase = read.delim("https://raw.githubusercontent.com/kkupkova/Mouse-bone-markers/main/BONE_DATABASE.tsv")
markerDatabase = markerDatabase[, !(names(markerDatabase) == "type")]
# encode greek letters to UTF-8
markerDatabase$gene = enc2utf8(markerDatabase$gene)
usethis::use_data(markerDatabase, overwrite = TRUE)
