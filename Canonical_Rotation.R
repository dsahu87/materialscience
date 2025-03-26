#Before use the script install "svDialogs" packages of R
#Vertical CO2 rotation
# Clean the strore data for a new start
rm(list = ls())
# Put the coordinates in the dialbox as it will be asked. 
{
  Cxx <-readline(prompt="coordinate Cx? ")
  Cyy <-readline(prompt="coordinate Cy? ")
  Czz <-readline(prompt="coordinate Cz? ")
  O1xx <-readline(prompt="coordinate O1x? ")
  O1yy <-readline(prompt="coordinate O1y? ")
  O1zz <-readline(prompt="coordinate O1z? ")
  O2xx <-readline(prompt="coordinate O2x? ")
  O2yy <-readline(prompt="coordinate O2y? ")
  O2zz <-readline(prompt="coordinate O2z? ")
Cx <- as.numeric(Cxx)
Cy <- as.numeric(Cyy)
Cz <- as.numeric(Czz)
O1x <- as.numeric(O1xx)
O1y <- as.numeric(O1yy)
O1z <- as.numeric(O1zz)
O2x <- as.numeric(O2xx)
O2y <- as.numeric(O2yy)
O2z <- as.numeric(O2zz)
## Relative coordinates with respect to Oxygen atom (Especially for cone like rotation)
N1x <- Cx-O2x
N1y <- Cy-O2y
N1z <- Cz-O2z
N2x <- O1x-O2x
N2y <- O1y-O2y
N2z <- O1z-O2z
# Mention angle in degree unit
A <- readline(prompt="rotation angle (A) ?" )  # angle of rotation in degree
#Convert all the angles in radian
split_num <-strsplit(A, ",")
angle<-as.numeric(c(unlist(split_num)))
Rad<-angle * pi / 180
#a <- as.numeric (A) * pi / 180 # convert angle in radians
x <- N1x 
y <- N1y 
z <- N1z
m <- N2x
n <- N2y
o <- N2z
##Rotation along the X-axis

#Getting relative Y1 cartesian coordinate after rotation for 1st rotating atom)
mult_1<-function(Rad,y,z)
{
  y*cos(Rad)-z*sin(Rad)
}
Y1<- mapply(mult_1,Rad,y,z)
#Getting relative Z1 cartesian coordinate after rotation for 1st rotating atom)
mult_2<-function(Rad,y,z)
{
  z*cos(Rad)+y*sin(Rad)
}
Z1<- mapply(mult_2,Rad,y,z) 
#Getting relative Y2 cartesian coordinate after rotation for 2nd rotating atom)
mult_3<-function(Rad,n,o)
{
  n*cos(Rad)-o*sin(Rad)
}
Y2<- mapply(mult_3,Rad,n,o)
#Getting relative Z2 cartesian coordinate after rotation for 2nd rotating atom)
mult_4<-function(Rad,o,n)
{
  o*cos(Rad)+n*sin(Rad)
}
Z2<- mapply(mult_4,Rad,o,n)

# C & O1 new coordinate after rotation [Oxygen O2 is freezed in vertical CO2]
NCy <-Y1+O2y
NCz <-Z1+O2z
NO1y <-Y2+O2y
NO1z <-Z2+O2z

data.frame(angle,NCy,NCz,NO1y,NO1z)
}
