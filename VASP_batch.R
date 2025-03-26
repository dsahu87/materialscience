## R code for performing DFT calculations on batch POSCAR files
## By Daniel Packwood (29-1-19)

## All POSCAR files must have the same elements in the same order.
## Identical KPOINTS, INCAR, and POTCAR files are used.

## OUTCAR files are analyzed and energies are automatically output

## POSCAR files must be named as POSCAR_k.vasp (k is the number)

rm(list = ls())

## 1. Basic input and output

input.folder <- "" # Location of the POSCAR files

calculation.folder <- paste(input.folder, "Calculations/", sep = "") # Directory where the calculation is performed. This folder must contain INCAR, KPOINTS, and POTCAR files.

output.folder <- paste(input.folder, "", sep = "") # Directory where OUTCAR files will be copied

dir.create(output.folder, recursive = TRUE, showWarnings = FALSE)

autoRunDFT = TRUE # Set to TRUE if DFT calculations are to be ran automatically

## 2. Sort the POSCAR files

batch.names <- dir()#dir(input.folder) # Get the names of the POSCAR files
poscar.flag <- grepl("POSCAR_", batch.names)

batch.name.ind <- c() # Get the indices of the folder names
for(k in batch.names[poscar.flag])
	{indk <- strsplit(k,"_c")[[1]][2];
	 indkj <- strsplit(indk, ".vasp")[[1]]

	 batch.name.ind <- c(batch.name.ind, as.numeric(indkj))
	}

batch.order <- sort(batch.name.ind,index.return=TRUE)$ix

poscar.names <- batch.names[poscar.flag][batch.order] # Puts the batch into correct order

## 3. Functions for the DFT calculation

base.dir <- getwd()

do.dft <- function(k) # Perform the DFT calculation for poscar.names[k]
	{command1 <- paste("cp ", input.folder, poscar.names[k], " ", calculation.folder, "POSCAR", sep = "")
	 system(command1)

	 setwd(calculation.folder)
		system("mpirun -np 32 vasp_std > log", wait = TRUE)
	 setwd(base.dir)

	 command2 <- paste("cp ", calculation.folder, "OUTCAR", " " , output.folder, "OUTCAR_", k, sep = "")
	 system(command2, wait = TRUE)	

	 command3 <- paste("cp ", calculation.folder, "CONTCAR", " " , output.folder, "CONTCAR_", k, sep = "")
	 system(command3, wait = TRUE)
	}

if(autoRunDFT)
	{for(k in 1:length(poscar.names))
		{do.dft(k)
		 print(paste("Job ", k, " of ", length(poscar.names), " is complete", sep = ""))
		}
	}

