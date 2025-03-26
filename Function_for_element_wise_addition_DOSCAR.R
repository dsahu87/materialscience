

## To use, edit the string "file.directory" as appropriate (see section 1 below), then run entire script in R or R Studio. Type pdos(k) to plot partial density of states for element k. Black line = total DOS, red line = s-orbital partial DOS, green line = p-orbital partial DOS, blue line = d-orbital partial DOS.

rm(list = ls())

## 1. Basic input (this part should be modified)


CONTCAR <- file.choose()
DOSCAR <-file.choose()

## 2. Read element information from POSCAR file (parts from here should not be modified)

CONT <- readLines(CONTCAR)
DOS <- readLines(DOSCAR)

elements <- regmatches(CONT[6], gregexpr("[^ ]+", CONT[6])) # Types of elements in the POSCAR file
n.elements <- regmatches(CONT[7], gregexpr("[^ ]+",CONT[7])) # Numbers of atoms for each element
n<-as.numeric(unlist(n.elements))
## 3. Read total DOSCAR information
dos.header<- regmatches(DOS[6], gregexpr("[^ ]+", DOS[6]))
head<-as.numeric(unlist(dos.header))
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

block<-function(k)
{
  ch<-0
  if (k>1)
    ch<-(sum(n[1:(k-1)])*ndiv)
  
  ch
  
  dat1<-table1[-(1:ch),]
  dat2<-dat1[1:(n[k]*ndiv),]
  pic<-as.data.frame(split(dat2, 1:n[k]))
  s.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V2", colnames(pic) )]) #extract all V2
s.orb}
  p1.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V3", colnames(pic) )]) #extract all V3
  p2.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V4", colnames(pic) )]) #extract all V4
  p3.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V5", colnames(pic) )]) #extract all V5
  d1.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V6", colnames(pic) )]) #extract all V6
  d2.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V7", colnames(pic) )]) #extract all V7
  d3.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V8", colnames(pic) )]) #extract all V8
  d4.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V9", colnames(pic) )]) #extract all V9
  d5.orb<-rowSums(pic[grepl( "^X\\d{1,3}.\\V10", colnames(pic) )]) #extract all V10
  
  plot(eax,dos,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)))
  par(new=TRUE) 
  plot(eax,s.orb,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=3)
  par(new=TRUE)
  plot(eax,(p1.orb+p2.orb+p3.orb),type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=2)
  par(new=TRUE)
  plot(eax,(d1.orb+d2.orb+d3.orb+d4.orb+d5.orb),type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)), col=4)
  abline(v = 0, lty = 2)
}
# Elementwise splitting data
#page(tot) to see full data
tot<-table1[1:57792,]
pic<-as.data.frame(split(tot, 1:192))
# pic[grepl( "^X\\d{1,3}.\\V1", colnames(pic) )] #extract all V1
s.orb<-pic[grepl( "^X\\d{1,3}.V2", colnames(pic) )] #extract all V2
p1.orb<-pic[grepl( "^X\\d{1,3}.V3", colnames(pic) )] #extract all V3
p2.orb<-pic[grepl( "^X\\d{1,3}.V4", colnames(pic) )] #extract all V4
p3.orb<-pic[grepl( "^X\\d{1,3}.V5", colnames(pic) )] #extract all V5
d1.orb<-pic[grepl( "^X\\d{1,3}.V6", colnames(pic) )] #extract all V6
d2.orb<-pic[grepl( "^X\\d{1,3}.V7", colnames(pic) )] #extract all V7
d3.orb<-pic[grepl( "^X\\d{1,3}.V8", colnames(pic) )] #extract all V8
d3.orb<-pic[grepl( "^X\\d{1,3}.V9", colnames(pic) )] #extract all V9
d3.orb<-pic[grepl( "^X\\d{1,3}.V10", colnames(pic) )] #extract all V10

lv<-split(tot,rep(1:192,each=ndiv))
gg<-as.data.frame(lv)
#rowSums(gg[grepl("^X\\d{1,3}", colnames(gg))])
s.orb<-rowSums(gg[grepl("^X\\d{1,3}.\\V2", colnames(gg))])

rowSums[grepl( "^X\\d{1,3}.V3", colnames(pic) )] # sum all V3 in row wise