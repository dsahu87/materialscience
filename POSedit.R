 # set the working directory
setwd("/Users/dsahu/Desktop/")

# "file.choose()" can be used here to choice entry
file <- "/Users/dsahu/Desktop/POSCAR" 

# read the file (header lines 9 are keeping hidden)
filedata <-read.table(file, skip=9, colClasses=c("character"))

# which data want to change and how much in detail
changedata <- as.numeric(filedata[192:195,"V3"]) + 1

# paste the changes data to the specific positions
filedata[192:195,"V3"]<- paste(changedata)

# save the file
write.table(filedata, file="/Users/dsahu/Desktop/POSCAR_modified")
