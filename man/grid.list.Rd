The term  grid.list refers to a list that contains information for 
evaluating a function on a 2-dimensional grid of points. If a function 
has more than two independent variables then one also needs to specify 
the constant levels for the variables that are not being varied. Also 
one needs to be able to  specify which variable should be identified 
with the horizontal or "X"  axis and which variable should be associated 
with the vertical or "Y" axis when the surface is 
plotted using contour or persp.  


This format is used extensively in the FUNFITS 
functions that evaluate function estimates and plot surfaces.


The form of a grid.list is 


list( var.name1= what  , var.name2=what , \dots{} var.nameN=what) 


Here var.names are the names of the independent variables.


The what options describe what should be done with this variable when 
generating the grid. The options are numerous reflecting some laziness 
in not wanting to type too much.  


'what'  can take on the following values:


"x" Generate the grid using this variable as the  horizontal axis  
for  a contour plot.


"y" Take this variable to be the vertical axis for a contour plot.


Numerical vector with two or more elements.  If only two elements these 
are interpreted as the _range_ for the grid. ( make.surface.grid has an 
option for specifying the number of grid points in this range with the 
default being 30.) If the vector has more than two elements these are 
assumed to be the actual grid values. Most of time  these are equally 
spaced and are easily generated using the seq function.  If two vectors 
are used to specify the x and y grids the first variable mentioned in 
the list is assumed as the x coordinate for contour and the second 
variable in order is the y coordinate. (There is an option in the 
as.surface function to flip these before plotting.) 


"c" This variable is held constant at its mean level. 


Numerical value.  This variable held constant at the specified value.


There are are some important defaults that are assumed in this format.
If a variable is omitted in the grid.list it is assumed to have value 
"c". If the variable names are omitted then the components are assumed 
to be in the same order as the variables. In this case one is not 
legal to omit a variable. Finally if a plotting or predict 
function takes a grid.list as an 
argument the default value for when it is not passed is list("x","y").
\seealso{
make.surface.grid, predict.surface, plot.surface, surface
}
\examples{
Given below are some examples of grid.list and the results 
when they are used with make.surface.grid. Note that 
make.surface.grid returns a matrix that retains the grid.list 
information as an attribute.


grid.l<- list( 1:3, 2:3)
make.surface.grid(grid.l) 
grid.l <- list( 1:3, 10, 1:3)
make.surface.grid(grid.l)


The next set of examples show how the grid.list can be used to 
control surface plotting and evaluation of an estimated function.


# first create a test function 
X<- 2*cbind( runif(50), runif(50), runif(50)) 
dimnames( X)<- list(NULL, c("X1","X2","X3"))
y<- X[,1]**2 + X[,2]**2 + exp(X[,3])  
# fit an  interpolating thin plate spline 
out<- tps( X,y,0)


grid.l<- list( X1=x, X2=y)
surface( out, grid.l) # surface plot of estimated surface with 
#                       X3= mean value


grid.l<- list( X1= seq( 0,1,,25), X2=0.0, X3=seq(0,1,,25))
surface( out,grid.l) #  surface plot based on a 25X25 grid in X1 an X3 
#                       over the square [0,1] and [0,1]  
#                       holding X2 equal to 0.0.
}
\keyword{FUNFITS}
% Converted by Sd2Rd version 0.2-a3.






          








