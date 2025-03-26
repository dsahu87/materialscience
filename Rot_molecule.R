#Before use the script install "svDialogs" packages of R
# Horizontal CO2 rotation
# Clean the strore data for a new start
rm(list = ls())
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
# Relative coordinates with respect to Carbon atom 
N1x <- O1x-Cx
N1y <- O1y-Cy
N1z <- O1z-Cz
N2x <- O2x-Cx
N2y <- O2y-Cy
N2z <- O2z-Cz
# Mention angle in degree unit
A <- readline(prompt="rotation angle (A)? ")  # angle of rotation in degree
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
##Rotation along the Z-axis

#Getting X1 values
mult_1<-function(Rad,x,y)
{
  x*cos(Rad)-y*sin(Rad)
}
X1<- mapply(mult_1,Rad,x,y)
#Getting Y1 values
mult_2<-function(Rad,x,y)
{
  y*cos(Rad)+x*sin(Rad)
}
Y1<- mapply(mult_2,Rad,x,y)
#Getting X2 values
mult_3<-function(Rad,m,n)
{
  m*cos(Rad)-n*sin(Rad)
}
X2<- mapply(mult_3,Rad,m,n)
#Getting Y2 values
mult_4<-function(Rad,m,n)
{
  n*cos(Rad)+m*sin(Rad)
}
Y2<- mapply(mult_4,Rad,m,n)

# Oxygen atoms (O1 & O2) new cartesian coordinate after rotation[carbon is freezed in horizontal CO2]
NO1x <- X1+Cx
NO1y <- Y1+Cy
NO2x <- X2+Cx
NO2y <-Y2+Cy
# Show the new coordinate values in a tabulated way
data.frame(angle,NO1x,NO1y,NO2x,NO2y)
}