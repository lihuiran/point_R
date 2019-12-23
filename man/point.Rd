\name{findpoints}
\alias{findpoints}
\title{label}
\usage{
find the label of the points which are classify as tip, branch and intermediate
}
\description{
return the label of each points
}
\examples{
library(matlab)

X = read.csv('C:/Users/dell/Desktop/Rfindpoints/data.csv',header=FALSE)

a = label(X,20,3,0.26)
}
