# R script for plotting density of states from DOSCAR file
# Start: 6-9-18

rm(list = ls())

# 1. Open DOSCAR file

dos.file <- "/Users/dsahu/Desktop/DOS/CO2_Fe-Pn_on_Cu_3_layer/DOSCAR" # Location of the DOSCAR file
nAtom <- 8 # Number of atoms in the unit cell

DOSdata <- readLines(dos.file) # Open the DOSCAR file

# 2. Function to extract numbers from a string 'line'

get.numbers <- function(line)
	{linesplit <- strsplit(line, " ")[[1]]
	 
	 linesplit <- as.numeric(linesplit)

	 whN <- which(!is.na(linesplit))

	 linesplit[whN]
	}

# 3. Get header information

header <- get.numbers(DOSdata[6])

Emin <- header[2] # Minumum energy in the DOS range
Emax <- header[1] # Maximum energy in the DOS range
nDiv <- header[3] # Number of energy divisions across the DOS range
Efermi <- header[4] # Fermi energy

# 4. Extract the DOS data

DOSlines <- DOSdata[(6+1):(6+nDiv)]

Ene <- c()
DOS <- c()

for(k in 1:nDiv)
	{DOSk <- get.numbers(DOSlines[k])

	 Ene <- c(Ene, DOSk[1])
	 DOS <- c(DOS, DOSk[2])
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(-15,Emax),ylim=c(0,100)) # Plot the density of states

abline(v = Efermi, lty=3) # Plot a vertical line at the Fermi level

# plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(-15,Emax),ylim=c(0,2000))
abline(v = Efermi, lty=3)

# 5. Plot the projected DOS data for the 1st atom

pDOSlines <- DOSdata[(6+nDiv+1+1):(6+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35)) # Plot the total density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=4) # Plot the d-contribution to the density of states

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),col=4) # Plot the d-contribution to the

# 6. Plot the projected DOS data for the 2nd atom.
pDOSlines <- DOSdata[(6+nDiv+1+nDiv+1+1):(6+nDiv+1+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35)) # Plot the total density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=4) # Plot the d-contribution to the density of states

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),axes=FALSE,col=4) # Plot the d-contribution to the density of states

# 7. Plot the projected DOS data for the 3rd atom.
pDOSlines <- DOSdata[(6+nDiv+1+nDiv+1+nDiv+1+1):(6+nDiv+1+nDiv+1+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35)) # Plot the total density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=4) # Plot the d-contribution to the density of states

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=4) # Plot the d-contribution to the density of states

# 8. Plot the projected DOS data for the 4th atom.
pDOSlines <- DOSdata[(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+1):(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35)) # Plot the total density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=4) # Plot the d-contribution to the density of states

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,15),axes=FALSE,col=4) # Plot the d-contributi

# 9. Plot the projected DOS data for the 5th atom.
pDOSlines <- DOSdata[(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+1):(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35)) # Plot the total density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,35),axes=FALSE,col=4) # Plot the d-contribution to the density of states

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=4) # Plot the d-contributi

# 10. Plot the projected DOS data for the 6th atom.
pDOSlines <- DOSdata[(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+1):(6+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv+1+nDiv)]

sDOS <- c() # s-orbital contribution to the DOS
pDOS <- c() # p-orbital contribution to the DOS
dDOS <- c() # d-orbital contribution to the DOS

for(k in 1:nDiv)
	{pDOSk <- get.numbers(pDOSlines[k])
	 
	 s <- pDOSk[2] # s-orbital contribution
	 p <- sum(pDOSk[3:5]) # Sum of p-orbital contributions
	 d <- sum(pDOSk[6:10]) # Sum of d-orbital contributions
	 
	 sDOS <- c(sDOS, s)
	 pDOS <- c(pDOS, p)
	 dDOS <- c(dDOS, d)
	}

plot(Ene,DOS,type='h',xlab="Energy (eV)", ylab="Number of states per unit cell",lwd=2,xlim=c(Emin,Emax),ylim=c(0,1300)) # Plot the total density of states
#par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
#plot(Ene,sDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,1300),axes=FALSE,col=2) # Plot the s-contribution to the density of states
#par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
#plot(Ene,pDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,1300),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS*100,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(Emin,Emax),ylim=c(0,1300),axes=FALSE,col=4) # Plot the d-contribution to the density of states

abline(v = Efermi, lty=3) # Plot a vertical line at the Fermi level

#plot(Ene,sDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),col=2) # Plot the s-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,pDOS*10,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=3) # Plot the p-contribution to the density of states
par(new = TRUE) # Allow for a new plot to be written on top of exisiting plot
plot(Ene,dDOS,type='h',xlab=" ", ylab=" ",lwd=2,xlim=c(-15,Emax),ylim=c(0,10),axes=FALSE,col=4) # Plot the d-contributi


