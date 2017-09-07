# R script tests assignPointsToExistingClusters(...) (APTEC) on synthesized data.

# Created by Sean Hendryx
# seanmhendryx@email.arizona.edu https://github.com/SMHendryx/assignPointsToClusters
# Copyright (c)  2017 Sean Hendryx
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################################################################################################

library(data.table)
library(ggplot2)
#library(truncnorm)
#library(tmvtnorm)


generatePointsInCircle = function(numSamples, min = 0, max = 1, plotter = FALSE){
  # Generates numSamples random points from a uniform distribution, within the radius between [min, max]
  # Returns a dataframe of x, y points 
  # adapted from: https://stackoverflow.com/questions/8473056/how-do-i-plot-a-circle-with-points-inside-it-in-r
  r = runif(numSamples, min = min, max = max)
  degs = 360*runif(numSamples)

  # First you want to convert the degrees to radians
  theta = 2*pi*degs/360

  # Plot your points by converting to cartesian
  xs = r*sin(theta)
  ys = r*cos(theta)
  X = cbind(xs, ys)
  colnames(X) = c('x', 'y')

  #plot:
  if(plotter){
    plot(X, xlim=c(-max(r),max(r)),ylim=c(-max(r),max(r)), asp = TRUE, main = "Generated Points")
    # Add a circle around the points
    polygon(max(r)*sin(seq(0,2*pi,length.out=100)),max(r)*cos(seq(0,2*pi,length.out=100)))
  }

  return(X)
}

generatePointOnCircle = function(center,r){
  # generates a point on the circle with a center at the two element vector center with radius r.
  numPoints = 1
  
  angle = runif(numPoints,0)

  x = center[1] + (r * cos(angle))
  y = center[2] + (r * sin(angle))

  point = c(x,y)

  return(point)
}

#test generatePointOnCircle():
r = 1
testCenter = c(2,2)
testPoint = generatePointOnCircle(center, r)
sqrt((testPoint[1]-testCenter[1])^2 + (testPoint[2]-testCenter[2])^2) == r


t1 = generatePointsInCircle(numSamplesPerTree, plotter = TRUE) + t1Centroid
dev.new()
plot(t1, asp = TRUE)

#generate clusters with true cluster identity (cluster) and a made-up clustering (clusterHat):
t1Centroid = c(2,2)
t2Centroid = c(4,2)
t3Centroid = c(7,4)

#generate sample points using truncated norm:
numSamplesPerTree = 1000
radius = 1.0

dt = data.table(x = numeric(), y = numeric(), tree = character())
trees = c('t1', 't2', 't3')
centroids = c(c(2,2), c(4,2), c(7,4))
centroids = matrix(centroids, nrow = length(trees), byrow = TRUE)
i = 1
for(tree in trees){
  rm(X)
  print(tree)
  X = generatePointsInCircle(numSamplesPerTree)
  X = t(t(X) + centroids[i,])
  X = as.data.table(X)
  X$tree = tree
  dt = rbind(dt, X)
  i = i + 1
}

# Now generate the points that are to be assigned to existing clusters:
points = data.table(x = rep(0.0, length(trees)), y = rep(0.0, length(trees)), tree = trees)
r = 1.0
i = 1
setkey(points, tree)
for(tree in trees){
  point = generatePointOnCircle(centroids[i,], r = r)
  points[tree, x := point[1]]
  points[tree, y := point[2]]
  i = i + 1
}

# add noise
noisyPoints = copy(points)
noisyPoints[,x := x + rnorm(3, sd = .05)]
noisyPoints[,y := y + rnorm(3, sd = .05)]

# plot:
p = ggplot(data = dt, mapping = aes(x, y, color = tree)) + geom_point() + theme_bw() + coord_equal()

# Add cluster labels
dt[, Label := ifelse(tree == 't3', 'c2', 'c1')]

p2 = ggplot(data = dt, mapping = aes(x, y, color = Label)) + geom_point() + theme_bw() + coord_equal()
#add noisyPoints:
p2 = p2 + geom_point(data = noisyPoints, mapping = aes(x = x, y = y, color = tree), shape = 8, size = 5)
p2 = p2 + guides(colour = guide_legend(override.aes = list(shape = c(16,16,8,8,8))))

# Now structure data to run through APTEC:

# Run APTEC on synthesized points and clusters:





#makes square:
genSamples_rtruncnorm = function(numSamplesPerTree, radius){
  t1x = t1Centroid[1] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t1y = t1Centroid[2] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t1 = cbind(t1x, t1y)
  colnames(t1) = c('x', 'y')
  t1 = as.data.table(t1)
  t1$tree = 't1'

  t2x = t2Centroid[1] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t2y = t2Centroid[2] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t2 = cbind(t2x, t2y)
  colnames(t2) = c('x', 'y')
  t2 = as.data.table(t2)
  t2$tree = 't2'

  t3x = t3Centroid[1] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t3y = t3Centroid[2] + rtruncnorm(numSamplesPerTree, a=-radius, b=radius, mean = 0, sd = 1)
  t3 = cbind(t3x, t3y)
  colnames(t3) = c('x', 'y')
  t3 = as.data.table(t3)
  t3$tree = 't3'

  dt = rbindlist(list(t1, t2,t3))

  p = ggplot(data = dt, mapping = aes(x, y, color = tree)) + geom_point() + theme_bw() + coord_equal()
}
