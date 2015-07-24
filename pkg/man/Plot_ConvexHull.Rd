\name{Plot_ConvexHull}
\alias{BMhyd}
\title{
	Convex hull
	}
	
\description{
	Plotting the convex hall
	}

\usage{
	Plot_ConvexHull(xcoord, ycoord, lcolor)
	\arguments{
  		\item{xcord}{coordinate for the horizontal axis}
  		\item{ycord}{coordinate for the vertical axis}
  		\item{col}{color}
  		}

\details{
	The Plot_ConvexHall function utilize the base function 'chull()', which takes x,y coordinates of your data and outputs a vector(hpts variable) of integers that indicate the points in your dataset that comprise the convex hull around your data. Then uses the base function lines to draw the convex hull using teh x,y coordinates stored in hpts.
	} 

\value{
	The convex hull plot for a set of 2D points.
 	}

\author{
	Ken Takagi 
	}

\references{
	https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/ 
	}

\examples{
	# Create a set of random data to plot convex hull around
	x<-rnorm(100,0.8,0.3)
	y<- rnorm(100,0.8,0.3)
	#get max and min of all x and y data for nice plotting
	xrange<-range(x)
	yrange<-range(y)
	#plot it up!
	plot(x,y,type="p",pch=1,col='black',xlim=c(xrange),ylim=c(yrange))
	Plot_ConvexHull(xcoord=x,ycoord=y,lcolor='black')
	}