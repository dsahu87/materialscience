setwd("/Users/dsahu/Desktop/")
file <- "/Users/dsahu/Desktop/POSCAR"
filedata <-read.table(file, skip=9, colClasses=c("character")) 
changedata <- as.numeric(filedata[192:195,"V3"]) + 1
filedata[192:195,"V3"]<- paste(changedata)
write.table(filedata, file="/Users/dsahu/Desktop/POSCAR_modified")