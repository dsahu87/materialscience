## R script for reading partial density-of-states, output from VASP.
## By Daniel Packwood (start: 8-1-19)

## POSCAR and DOSCAR must be located in the same folder!
## POSCAR must be in the format

## CuSCN (Smith 1982)                      
##   1.00000000000000     
##     3.8621779709993063   -0.0000000000000000   -0.0000000000000000
##    -1.9310889854996531    3.3447442368712355    0.0000000000000000
##     0.0000000000000000    0.0000000000000000   10.9670315265920060
##   Cu   N    C    S 
##     2     2     2     2
## (etc)

## To use, edit the string "file.directory" as appropriate (see section 1 below), then run entire script in R or R Studio.
## Type pdos(k) to plot partial density of states for element k. 
## Black line = total DOS, red line = s-orbital partial DOS, green line = p-orbital partial DOS, blue line = d-orbital partial DOS.

rm(list = ls())

## 1. Basic input (this part should be modified)

file.directory <- "/Users/dsahu/Desktop/DOS/Double_Pn/T3_VER_CO2_with_OPT_Cu_12_5_1_added2Pn_FAR_fromOUTPUT/" # Directory where POSCAR and DOSCAR files are located

## 2. Read element information from POSCAR file (parts from here should not be modified)

POSCAR <- readLines(paste(file.directory,"POSCAR",sep=""))

get.elements <- function(POSCAR) # Extract element types from POSCAR
{eline <- strsplit(POSCAR[6], " ")[[1]]

eline[eline!=""]
}

get.numbers <- function(line) # Extract numbers from the string 'line'
{nline <- strsplit(line, " ")[[1]]
nline <- as.numeric(nline)

nline[!is.na(nline)]
}

elements <- get.elements(POSCAR) # Types of elements in the POSCAR file
n.elements <- get.numbers(POSCAR[7]) # Numbers of atoms for each element

## 3. Read total DOSCAR information

DOSCAR <- readLines(paste(file.directory, "DOSCAR", sep = ""))

dos.header <- function(DOSCAR) # Extract the header information (energy range, number of divisions, Fermi energy) from the DOSCAR file
{head <- get.numbers(DOSCAR[6])

c(head[2],head[1],head[3],head[4])
}

dhead <- dos.header(DOSCAR)

emin <- dhead[1]
emax <- dhead[2]
ndiv <- dhead[3]
efermi <- dhead[4]

get.dos <- function(DOSCAR) # Extract the total DOS from the DOSCAR file
{indS <- 7
indF <- 6 + ndiv

dat.out <- sapply(DOSCAR[indS:indF], get.numbers, USE.NAMES=FALSE)

list(dat.out[1,],dat.out[2,])
}

dost <- get.dos(DOSCAR)

eax <- dost[[1]] - efermi # Energy axis for plotting density of states. 0 will be the Fermi level
dos <- dost[[2]] # Total density of states

## 4. Functions to read and plot partial DOS for a given element type

get.index.at <- function(i) # Get the lines of DOSCAR for atom i
{i.start <- 6 + (ndiv + 1) + (ndiv + 1)*(i - 1) + 1
i.end <- 6 + (ndiv + 1) + (ndiv + 1)*(i - 1) + ndiv

i.start:i.end
}

get.index.el <- function(k) # Get the lines of DOSCAR for all atoms of type k
{atn.s <- 1
atn.f <- n.elements[1]

if(k > 1)
{atn.s <- sum(n.elements[1:(k-1)]) + 1
atn.f <- sum(n.elements[1:k])
}

atn <- atn.s:atn.f;

dos.index <- sapply(atn, get.index.at, USE.NAMES = FALSE) # Columns are indices for each atom

c(dos.index) 
}

pdos <- function(k) # Plot the partial density of states for element k
{
  pdos.s <- rep(0,ndiv) # p-dos for s-states
  pdos.p <- rep(0,ndiv) # p-dos for p-states
  pdos.d <- rep(0,ndiv) # p-dos for d-states
  
  dos.index <- get.index.el(k)
  
  pdos.dat <- sapply(DOSCAR[dos.index], get.numbers, USE.NAMES = FALSE)
  
  index.at <- list() # Assign the columns of pdos.dat to atoms
  count <- 1
  
  for(i in 1:n.elements[k])
  {iadd <- count:((i-1)*ndiv + ndiv)
  index.at[[i]] <- iadd
  
  count <- count + ndiv
  }
  
  # Sum the p-dos over the atoms
  
  for(i in 1:length(index.at))
  {pdos.s <- pdos.s + pdos.dat[2,index.at[[i]]]
  
  pdos.p <- pdos.p + pdos.dat[3,index.at[[i]]] + pdos.dat[4,index.at[[i]]] + pdos.dat[5,index.at[[i]]]
  
  pdos.d <- pdos.d + pdos.dat[6,index.at[[i]]] + pdos.dat[7,index.at[[i]]] + pdos.dat[8,index.at[[i]]] + pdos.dat[9,index.at[[i]]] + pdos.dat[10,index.at[[i]]]
  }
 
  pDOSlines <- DOSCAR[(6+1+((ndiv+1)*213)):(6+ndiv+((1+ndiv)*213))]
  sDOS <- c() # s-orbital contribution to the DOS
  pDOS <- c() # p-orbital contribution to the DOS
  dDOS <- c() # d-orbital contribution to the DOS
  
  for(k in 1:ndiv)
  {pDOSk <- get.numbers(pDOSlines[k])
  
  s <- pDOSk[2] # s-orbital contribution
  p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
  d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
  
  sDOS <- c(sDOS, s)
  pDOS <- c(pDOS, p)
  dDOS <- c(dDOS, d)
  }
  
  # Perform the plotting
  
  plot(eax,dos,type='l',lwd=2,xlab="energy (eV)",ylab="DOS",ylim=c(0,max(dos)))
  
  par(new=TRUE) 
  plot(eax,pdos.s,type='l',lwd=2,xlab="",ylab="",ylim=c(0,max(dos)),axes=FALSE,col=2)
  
  par(new=TRUE) 
  plot(eax,pdos.p,type='l',lwd=2,xlab="",ylab="",ylim=c(0,max(dos)),axes=FALSE,col=3)
  
  par(new=TRUE) 
  plot(eax,pdos.d,type='l',lwd=2,xlab="",ylab="",ylim=c(0,max(dos)),axes=FALSE,col=4)
  
  par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
  plot(eax,sDOS,type='l',lwd=2,xlab="", ylab="",xlim=c(emin,emax),ylim=c(0,max(dos)),axes=FALSE, col=2) # Plot the s-contribution to the density of states
  
  par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
  plot(eax,pDOS,type='l',lwd=2,xlab="", ylab="",xlim=c(emin,emax),ylim=c(0,max(dos)),axes=FALSE,col=3) # Plot the p-contribution to the density of states
  
  par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
  plot(eax,dDOS,type='l',lwd=2,xlab="", ylab="",xlim=c(emin,emax),ylim=c(0,max(dos)),axes=FALSE,col=4) # Plot the d-contribution to the density of states
  
  abline(v = 0, lty = 2)
  
  
}
