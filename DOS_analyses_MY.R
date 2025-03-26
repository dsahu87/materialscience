
## To use, edit the string "file.directory" as appropriate (see section 1 below), then run entire script in R or R Studio.
#Type block(k) to plot partial density of states for element k. Black line = total DOS, red line = s-orbital partial DOS, green line = p-orbital partial DOS, blue line = d-orbital partial DOS.

rm(list = ls())

## 1. Import CONTCAR & DOSCAR files as DATA

CONTCAR <- file.choose()
DOSCAR <-file.choose()

## 2. Read the imported DATA files CONTCAR & DOSCAR

CONT <- readLines(CONTCAR)
DOS <- readLines(DOSCAR)
## Extract the number of elements and their type
elements <- regmatches(CONT[6], gregexpr("[^ ]+", CONT[6])) # Types of elements in the POSCAR file
n.elements <- regmatches(CONT[7], gregexpr("[^ ]+",CONT[7])) # Numbers of atoms for each element
n<-as.numeric(unlist(n.elements))
## 3. Extracting the required information from DOSCAR
dos.header<- regmatches(DOS[6], gregexpr("[^ ]+", DOS[6])) # Extract data of line number 6 of DOSCAR with right way
head<-as.numeric(unlist(dos.header)) #Destroy the list style to treat as numeric
emax<-head[1]
emin<-head[2]
ndiv<-head[3]
efermi<-head[4]
skip<-DOS[-(1:308)]
skip1<-skip[-(grepl(paste(DOS[6]),skip,fixed=TRUE))==0]

table1<-read.table(textConnection(skip1),sep="")

hold<-DOS[(1:308)]
hold1 <- hold[-(1:grep("unknown system",hold))]
hold2<-hold1[-(grepl(paste(DOS[6]),hold1,fixed=TRUE))==0]
table2<-read.table(textConnection(hold2),sep="")
eax <- (table2$V1) - efermi # Energy axis for plotting density of states. 0 will be the Fermi level
dos <- (table2$V2) # Total density of states

## Single function to run all commands at once by typing block(k): k=1,2,3..etc as elemental order

block<-function(k)
  {
ch<-0
  if (k>1)
ch<-(sum(n[1:(k-1)])*ndiv) # Skip the number of lines before the 'K' element

ch

dat1<-table1[-(1:ch),]
dat2<-dat1[1:(n[k]*ndiv),] # skip the tail lines for the other elements 

lv<-split(dat2,rep(1:n[k],each=ndiv)) # split each elementwise for a element type
gg<-as.data.frame(lv) # Place all the splitted data in a frame
s.orb<-rowSums(gg[grepl("^X\\d{1,3}.\\V2", colnames(gg))]) # row wise sum of s-orb DOS for all elements with same type 
p1.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V3", colnames(gg))])# row wise sum of px-orb DOS for all elements with same type
p2.orb<-rowSums(gg[grepl("^X\\d{1,3}.\\V4", colnames(gg))])# row wise sum of py-orb DOS for all elements with same type
p3.orb<-rowSums(gg[grepl("^X\\d{1,3}.\\V5", colnames(gg))])# row wise sum of pz-orb DOS for all elements with same type
d1.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V6", colnames(gg))])# row wise sum of dxy-orb DOS for all elements with same type
d2.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V7", colnames(gg))])# row wise sum of dyz-orb DOS for all elements with same type
d3.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V8", colnames(gg))])# row wise sum of dxz-orb DOS for all elements with same type
d4.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V9", colnames(gg))])# row wise sum of dx2-y2-orb DOS for all elements with same type
d5.orb<- rowSums(gg[grepl("^X\\d{1,3}.\\V10", colnames(gg))])# row wise sum of dz2-orb DOS for all elements with same type
p.orb<-p1.orb + p2.orb + p3.orb # sum up px, py & pz orb DOS
d.orb<-d1.orb + d2.orb + d3.orb + d4.orb + d5.orb # sum up dxy, dyz, dxz, dx2-y2, dz2
plot(eax,dos,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos))) # Plot total DOS of full system
par(new=TRUE) 
plot(eax,s.orb*10,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=2) # Plot partial DOS for s-orb
par(new=TRUE)
plot(eax,p.orb*10,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=3)# Plot partial DOS for p-orb
par(new=TRUE)
plot(eax,d.orb*10,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=4)# Plot partial DOS for d-orb
abline(v = 0, lty = 2)
}


