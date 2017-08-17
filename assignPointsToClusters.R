# Algorithms for assigning points in one dataset to clusters in another dataset, flagging unlikely correspondences based on a distance threshold, 
# and "merging" clusters that have been over-segmented by the clustering algorithm.
# Ideally, if we have a dataset of points, where each point represents a cluster in another data set, there would be a one-to-one 
# correspondence between the two data sets, such that for each point there is one and only one cluster to which it corresponds and
# such that it is unmistakable which point corresponds to which cluster.  
# Though, this is not usually the case.
# 
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


#load packages:
library(data.table)


#helper functions first:
printer <- function(string, variable){
  print(paste0(string, variable))
}

#define Euclidean distance function for working with two points/vectors in n-dimensional space:
eucDist <- function(x1, x2){
  sqrt(sum((x1 - x2) ^ 2))
} 

#dmode function finds the mode with the highest peak (dominate mode) and nmodes identify the number of modes.
#https://stackoverflow.com/questions/16255622/peak-of-the-kernel-density-estimation
dmode <- function(x) {
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)]) 
}  

nmodes <- function(x) {  
  den <- density(x, kernel=c("gaussian"))
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.8)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  nmodes <- length(rle(den.sign <- sign(s.derv$s1))$values)/2
  if ((nmodes > 10) == TRUE) { nmodes <- 10 }
  if (is.na(nmodes) == TRUE) { nmodes <- 0 } 
  ( nmodes )
}

# Toy Example
#x <- runif(1000,0,100)
#plot(density(x))
#  abline(v=dmode(x))


# Algorithms:
####################################################################################################################################################################################

thresholdPoints <- function(points, thresholdType = "dominateMode", buffer = 10, plotDensity = FALSE){
  ######################################################################################################################################################
  # Function computes threshold beyond which it is unlikely that the point corresponds to the cluster.
  ######################################################################################################################################################
  # Add boolean column indicating if the closest cluster is beyond threshold
  #first copy points to be modified locally (not by reference)
  points = copy(points)
  points[, closest_cluster_outside_threshold := NULL]
  points[, closest_cluster_outside_threshold := logical()]
  # Return the distance of the closest in situ coordinate (i.e. "point" in points) to each cluster centroid in data.table:
  # .SD[] makes Subset of Datatable
  # https://stackoverflow.com/questions/33436647/group-by-and-select-min-date-with-data-table
  closestPoints = points[,.SD[which.min(distance_to_closest_cluster_member)], by = cluster_ID]

  # Compute the threshold beyond which it is unlikely that the point corresponds to the cluster 
  #   (or put another way, that the cluster represents the point):
  print(closestPoints)
  dominateModeDist = dmode(closestPoints$distance_to_closest_cluster_member)
  print(paste0("Dominate mode of distances between point and closest cluster member: ", dominateModeDist))
  meanDist = mean(closestPoints$distance_to_closest_cluster_member)
  print(paste0("Mean of distances between point and closest cluster member: ", meanDist))
  if (plotDensity) { #then:
    plot(density(closestPoints$distance_to_closest_cluster_member))
  }
  if(thresholdType == "dominateMode"){
    threshold = dominateModeDist * buffer
  }else{
    threshold = meanDist * buffer
  }
  print(paste0("threshold = ", threshold))
  points[, closest_cluster_outside_threshold := (distance_to_closest_cluster_member > threshold)]
  #
  return(points)
}

# vectorized:
assignPointsToExistingClusters <- assignPointsToClusters <- function(points, clusters, x_col_name = 'X', y_col_name = 'Y', cluster_ID_col_name = 'Label', thresholdType = "dominateMode", buffer = 10){
  # Algorithm assigns points, in a dataset $\bf{P}$, to the closet cluster in another dataset, $\bf{C}$,
  # flags unlikely correspondences based on distance threshold,
  # and then determines if any other clusters should be assigned to that point based on information held in the point.
    # if we have a matrix of point coordinates and each of the points represents a cluster, the algorithm assigns each point to a cluster.
  # if outliers are coded as -1 in cluster_ID_col_name, they will be assumed to not be clusters
  # :Param points: data.table object with columns 'X' and 'Y'
  # :Param clusters: data.table object with columns 'X', 'Y', and 'Label'
  # :Params x_col_name = 'X', y_col_name = 'Y', cluster_ID_col_name = 'Label': refer to columns in point data.table
  #check if column already exists:
  if(any(names(points) == "cluster_ID")){
    stop("cluster_ID already exists in points.  points should not include cluster ids prior to running assignPointsToClusters function.")
  }
  #first copy points and clusters to be modified locally (not by reference)
  points = copy(points)
  clusters = copy(clusters)
  #if doesn't exist, add:
  points[,cluster_ID := integer()]
  
  # for checkIfPointRepresentsMoreThanOneCluster
  points[,x_closestCentroid := double()]
  points[,y_closestCentroid  := double()]
  
  # remove outliers coded as -1:
  clusters = clusters[Label != -1,]
  clusterLabels = unique(clusters[,cluster_ID_col_name, with = FALSE])
  #print(paste0("clusterLabels", clusterLabels))

  # Now loop through points and find each point's closest cluster
  # because we are looping through the points and finding the closest cluster, multiple points can be assigned to the same cluster
  for(i in seq(nrow(points))){
    print(paste0("Looping through points to find closest cluster, on point: ", i))
    position = points[i, c(x_col_name, y_col_name), with = FALSE]
    #Trying to add position as a column:
    clusters[,X_pointPosition := position[[1]]]
    clusters[,Y_pointPosition := position[[2]]]
    #find shortest distance between point and any cluster member (cluster point):
    clusters[,XDiff := (X_pointPosition - X)]
    clusters[,YDiff := (Y_pointPosition - Y)]
    clusters[,Xsq := XDiff^2]
    clusters[,Ysq := YDiff^2]
    clusters[,summed := Xsq + Ysq]
    clusters[,distance_to_point := sqrt(summed)]
    closestMember = clusters[, .SD[which.min(distance_to_point)]]
    #
    print(paste0("closestMember: ", closestMember$Label))
    points[i,cluster_ID := closestMember$Label]
    print(paste0("closest cluster member distance = ", closestMember$distance_to_point))
    points[i,distance_to_closest_cluster_member := closestMember$distance_to_point]
    points[i, X_closest_cluster_member := closestMember$X]
    points[i, Y_closest_cluster_member := closestMember$Y]
    #compute the centroid of the closest cluster:
    X = clusters[Label == closestMember$Label, X]
    Y = clusters[Label == closestMember$Label, Y]
    clusterCentroid = colMeans(cbind(X, Y))
    points[i, X_closest_cluster_centroid := clusterCentroid[1]]
    points[i, Y_closest_cluster_centroid := clusterCentroid[2]]
  }
  points = thresholdPoints(points, thresholdType = thresholdType, buffer = buffer)
  return(points)
}

testIfPointWithinCircle <- function(x, center_x, y, center_y, radius){
  # tests if point falls within circle defined by a center point and a radius
  # ported from C, philcolbourn answer: https://stackoverflow.com/questions/481144/equation-for-testing-if-a-point-is-inside-a-circle
  dx = abs(x-center_x)
  dy = abs(y-center_y)
  #printer("dx: ", dx)
  if(dx > radius){
    return(FALSE)
  }
  if (dy > radius){
    return(FALSE)
  }
  if (dx + dy <= radius){
    return(TRUE)
  }
  if (dx^2 + dy^2 <= radius^2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


computeUnassignedClusterCentroids <- function(clusters){
  # Given a data.table of clusters, returns a data.table of the unassignedClusterCentroids
  # Make datatable of unassigned clusters:
  unassignedClusters = copy(clusters[is.na(assigned_to_point),])
  # Make copy of unassigned cluster centroids (to be filled in):
  unassignedClusterLabels = unique(unassignedClusters[,Label])
  # convert from factor to numeric:
  unassignedClusterLabels = as.numeric(unassignedClusterLabels)
  unassignedClusterCentroids = data.table(Label = unassignedClusterLabels)#, X =  double(), Y = double(), Z = double(), assigned_to_point = NA)
  setkey(unassignedClusterCentroids)
  unassignedClusterCentroids[,X:= double()]
  unassignedClusterCentroids[,Y:= double()]
  #just in 2D for now
  #unassignedClusterCentroids[,Z:= double()]
  unassignedClusterCentroids[,assigned_to_point := NA]
  
  # compute unassignedClusterCentroids:
  for (unassignedClusterLabel in unassignedClusterLabels){
    X = clusters[Label == unassignedClusterLabel, X]
    Y = clusters[Label == unassignedClusterLabel, Y]
    clusterCentroid = colMeans(cbind(X, Y))
    unassignedClusterCentroids[Label == unassignedClusterLabel, X := clusterCentroid[1]]
    unassignedClusterCentroids[Label == unassignedClusterLabel, Y := clusterCentroid[2]]
    #print("One iteration of for loop computing and storing unassigned cluster centroid.")
  }
  return(unassignedClusterCentroids)
}


testAndMergeClustersRecursively <- function(predictedCentroid, pointID, assignedPoints, clusters){
  # Function resolves over-segmentation of clusters by determining if any other clusters should be assigned to the points in assignedPoints
  # Updates clusters by reference (does not return a new object)
  # predictedCentroid: list of x, y predicted coordinate of center of true cluster (updated recursively)
  # pointID: string refers to the id of the point in assignedPoints for which we are searching for clusters that belong to it.
  # assignedPoints: data.table of points with assignments to clusters (values in assignedPoints$cluster_ID)
  # clusters: data.table containing clustered points
  # note that assignedPoints and clusters are two different datasets that represent the same things in the real world but have some small difference in their representation of the real world objects
  
  print("Testing clusters.")
  center_x = predictedCentroid[1]
  center_y = predictedCentroid[2]
  #center_x = assignedPoints[Sample_ID == pointID, X_closest_cluster_centroid]
  #center_y = assignedPoints[Sample_ID == pointID, Y_closest_cluster_centroid]
  radius = assignedPoints[Sample_ID == pointID, Minor_Axis]
  
  if(is.na(radius)){
    stop("radius does not exist.  Make sure all points have a radius value.")
  }
  
  #These two lines of code could only be run once and the data stored to make run faster.  i.e. compute unassigned labels, compute unassigned centroids => store unassigned centroids => remove unassigned centroid if the centroid becomes assigned.
  # Compute remaining unassigned cluster labels:
  unassignedClusterLabels = unique(clusters[is.na(assigned_to_point), Label])
  # Compute unassignedClusterCentroids:
  unassignedClusterCentroids = computeUnassignedClusterCentroids(clusters)
  
  for (unassignedClusterLabel in unassignedClusterLabels){
    # test if point falls within minor axis circle from assigned cluster centroid:
    #x, center_x, y, center_y, radius
    #print(paste0("Testing unassigned cluster: ", unassignedClusterLabel))
    x = unassignedClusterCentroids[Label == unassignedClusterLabel, X]
    y = unassignedClusterCentroids[Label == unassignedClusterLabel, Y]

    #print("testIfPointWithinCircle: ")
    #printer("x: ", x)
    #printer("center_x: ", center_x)
    #printer("y: ", y)
    #printer("center_y: ", center_y)
    #printer("radius: ", radius)
    if (testIfPointWithinCircle(x = x, center_x = center_x, y = y, center_y = center_y, radius = radius)){
      print(paste0("Cluster centroid falls within radius.  Cluster Label: ", unassignedClusterLabel))
      #if unassigned cluster centroid within minor_axis radius of assigned centroid,
      # assign point to cluster:
      clusters[Label == unassignedClusterLabel, assigned_to_point := assignedPoints[Sample_ID == pointID, Sample_ID]]
      # add column to indicate if cluster has been merged with another cluster to represent some single point:
      clusters[Label == unassignedClusterLabel, merged := TRUE]
      assignedPoints[Sample_ID == pointID, merged := TRUE]
      printer("pointID that represents more than one cluster: ", pointID)
      #printer("correspondence_ID: ", assignedPoints[Sample_ID == pointID, correspondence_ID])
      
      # Add new row recording the point and cluster correspondence
      newCorrespondenceID = (max(assignedPoints[,correspondence_ID]) + 1)
      #printer("correspondence_ID on newly assigned cluster: ", newCorrespondenceID)
      
      newAssignedPointsRow = copy(assignedPoints[Sample_ID == pointID])
      #printer("newAssignedPointsRow: ", newAssignedPointsRow)
      newAssignedPointsRow[,correspondence_ID := newCorrespondenceID]
      #printer("newAssignedPointsRow with changed correspondence_ID: ", newAssignedPointsRow)
      assignedPoints = rbindlist(list(assignedPoints, newAssignedPointsRow))
      #printer("assignedPoints[correspondence_ID == newCorrespondenceID,] ", assignedPoints[correspondence_ID == newCorrespondenceID,])
      
      # Add previously unassigned cluster label to assignedPoints
      assignedPoints[correspondence_ID == newCorrespondenceID, cluster_ID := unassignedClusterLabel]
      assignedPoints[correspondence_ID == newCorrespondenceID, merged := TRUE]
      #printer("after merged := true, assignedPoints[correspondence_ID == newCorrespondenceID,] returns:", assignedPoints[correspondence_ID == newCorrespondenceID,])
      
      
      # compute newPredictedCentroid from clusters:
      X = clusters[assigned_to_point == pointID, X]
      Y = clusters[assigned_to_point == pointID, Y]
      newPredictedCentroid = colMeans(cbind(X, Y))
      print(paste0("New predicted centroid: ", newPredictedCentroid))
      
      #Recursive call:
      print("Starting recursive call to testAndMergeClustersRecursively:")
      #nextAssignedPoints = copy(assignedPoints)
      assignedPoints = testAndMergeClustersRecursively(newPredictedCentroid, pointID, assignedPoints, clusters)
    }# end if (testIfPointWithinCircle)
  }# end for unassignedClusterLabel in unassignedClusterLabels
  return(assignedPoints)
}

checkIfPointRepresentsMoreThanOneCluster <- function(assignedPoints, clusters){
  # Function resolves over-segmentation of clusters by determining if any other clusters should be assigned to the points in assignedPoints
  # based on information held in the point (metadata) and, if so,
  # assigns the cluster(s) to the point.
  # Assumes that metadata is at the same scale of clusters
  # Returns assignedPoints with new column correspondence_ID, which contains the id of a unique correspondence relationship between a point and a cluster.
  
  print("Resolving over-segmentation of clusters.  Checking all points to see if any point represents more than one cluster.")
  
  #First add correspondence column:
  for(i in seq(nrow(assignedPoints))){
    assignedPoints[i,correspondence_ID := i]
    #printer("correspondence_ID: ", i)
  }
  
  # Updates the clusters data.table object by reference.  
  # First, remove clusters outliers coded as -1:
  clusters = clusters[Label != -1,]
  # if Label is factor, Label -1 still exists as level, so:
  #clusters[,Label := droplevels(Label)]
  
  # convert assignedPoints$Sample_ID from factor (likely default) to character:
  assignedPoints[,Sample_ID := as.character(Sample_ID)]
  
  # Add assignedPoint ID (currently hardcoded as Sample_ID) to clusters:
  # to vectorize, do something like this: clusters[,assigned_to_point := ].  Otherwise:
  for(i in seq(nrow(assignedPoints))){
    # Adding assignedPoint ID to clusters:
    point = copy(assignedPoints[i,])
    clusters[Label == point$cluster_ID, assigned_to_point := point$Sample_ID]
  }
  
  clusters[,assigned_to_point := as.character(assigned_to_point)] 
  
  # Now loop through assignedPoints, to see if any unassigned cluster centroids fall within Minor_Axis radius from assigned cluster centroid:
  # Making list of those points that are uniquely assigned to a cluster, as it is unlikely that the cluster is over-segmented if more than one point has been assigned to the cluster:
  uniquelyAssignedPointIDs = vector(mode = "character")
  assignedClusterIDs = unique(assignedPoints$cluster_ID)
  
  for (id in assignedClusterIDs){
    if(nrow(assignedPoints[cluster_ID == id]) == 1) {
      uniquelyAssignedPointIDs_i = assignedPoints[cluster_ID == id, Sample_ID]
      # append two vectors together:
      uniquelyAssignedPointIDs = c(uniquelyAssignedPointIDs, uniquelyAssignedPointIDs_i)
      #printer("uniquelyAssignedPointIDs: ", uniquelyAssignedPointIDs)
      print(uniquelyAssignedPointIDs)
    } 
  }
  
  for(pointID in uniquelyAssignedPointIDs){
    print(paste0("Testing if any other clusters belong to point: ", pointID))
    # compute predictedCentroid of true cluster from clusters (i.e., compute the centroid of the points in cluster datatable that are assigned to point with pointID):
    X = clusters[assigned_to_point == pointID, X]
    Y = clusters[assigned_to_point == pointID, Y]
    predictedCentroid = colMeans(cbind(X, Y))
    
    assignedPoints = testAndMergeClustersRecursively(predictedCentroid = predictedCentroid, pointID = pointID, assignedPoints = assignedPoints, clusters = clusters)
  }
  #remove duplicate rows:
  assignedPoints = unique(assignedPoints)
  # instead of returning clusters, add list of cluster_IDs column to assignedPoints and return assignedPoints:
  return(assignedPoints)
}



buildClusterDict <- function(assignedPoints, evalOutsideThreshold = FALSE){
  # Returns a dictionary where each sample ID (point) represents one to many clusters
  # takes in a datatable that has been run through assignPointsToExistingClusters(...) and checkIfPointRepresentsMoreThanOneCluster(...)
  # The dictionary is just a list of lists
  # where the name of each entry is the Sample_ID (training and validation ID) and the entries are the IDs of the clusters (cluster_ID in assignedPoints and Label in clusters)
  sampleIDs = as.list(unique(assignedPoints[closest_cluster_outside_threshold == evalOutsideThreshold,Sample_ID]))
  clusterDict = vector(mode = "list", length= length(sampleIDs))
  names(clusterDict) = sampleIDs
  i = 1
  for(sampleID in sampleIDs){
    clusterIDs_i = assignedPoints[Sample_ID==sampleID, cluster_ID]
    clusterIDs_i = as.list(clusterIDs_i)
    clusterDict[[i]] = clusterIDs_i
    i = i + 1 
  }
  return(clusterDict)
}

assignMergedIDsToClusters <- function(clusterDict, clusters){
  # Function assigns a new ID of merged clusters to clusters datatable
  # Run after buildClusterDict(...)
  clusters = copy(clusters)
  #note that pointID and sampleID (in points) represent the same thing in this code base
  pointIDs = names(clusterDict)
  for(pointID in pointIDs){
    entry = clusterDict[pointID]
    clusters[Label %in% entry[[1]],mergedClusterID := pointID]
  }
  return(clusters)
}

