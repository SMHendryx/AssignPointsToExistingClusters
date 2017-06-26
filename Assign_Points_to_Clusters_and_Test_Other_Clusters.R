# Algorithms for assigning points in one dataset to clusters in another dataset 
# flagging unlikely correspondences based on a distance threshold.
# Authored by Sean Hendryx while working at the University of Arizona

#load packages:
library(data.table)

#helper functions first:
#define Euclidean distance function for working with two points/vectors in n-dim space:
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


#algos:
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
assignPointsToClusters <- function(points, clusters, x_col_name = 'X', y_col_name = 'Y', cluster_id_col_name = 'Label', thresholdType = "dominateMode", buffer = 10){
  # Algorithm assigns points, in a dataset $\bf{P}$, to the closet cluster in another dataset, $\bf{C}$,
  # flags unlikely correspondences based on distance threshold,
  # and then determines if any other clusters should be assigned to that point based on information held in the point.
    # if we have a matrix of point coordinates and each of the points represents a cluster, the algorithm assigns each point to a cluster.
  # if outliers are coded as -1 in cluster_id_col_name, they will be assumed to not be clusters
  # :Param points: data.table object with columns 'X' and 'Y'
  # :Param clusters: data.table object with columns 'X', 'Y', and 'Label'
  #check if column already exists:
  if(any(names(points) == "cluster_ID")){
    stop("cluster_ID already exists in points.  points should not include cluster ids prior to running assignPointsToClusters function.")
  }
  #first copy points and clusters to be modified locally (not by reference)
  points = copy(points)
  clusters = copy(clusters)
  #if doesn't exist, add:
  points[,cluster_ID := integer()]
  # are these columns necessary????????????????????????????????????????????????????:
  # YES!  for checkIfPointRepresentsMoreThanOneCluster
  points[,x_closestCentroid := double()]
  points[,y_closestCentroid  := double()]
  # ^^^^^^^^^^^^^^^^^Necessary ^???????????????????????????????????????????????????
  # remove outliers coded as -1:
  clusters = clusters[Label != -1,]
  clusterLabels = unique(clusters[,cluster_id_col_name, with = FALSE])
  #print(paste0("clusterLabels", clusterLabels))
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


checkIfPointRepresentsMoreThanOneCluster <- function(assignedPoints, clusters){
  # Function determines if any other clusters should be assigned to the point based on information held in the point (metadata) and, if so,
  # assigns the cluster(s) to the point.
  # Assumes that metadata is at the same scale of clusters
  # returns the clusters
  
  # First, remove clustersoutliers coded as -1:
  clusters = clusters[Label != -1,]
  # since factor, Label -1 still exists as level, so:
  clusters[,Label := droplevels(Label)]
  # Add assignedPoint ID (currently hardcoded as Sample_ID) to clusters:
  # to vectorize, do something like this: clusters[,assigned_to_point := ].  Otherwise:
  for(i in seq(nrow(points))){
    # Adding assignedPoint ID to clusters:
    point = copy(assignedPoints[i,])
    clusters[Label == point$cluster_ID, assigned_to_point := point$Sample_ID]
  }#end for loop
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
  # I am here
  # compute unassignedClusterCentroids:
  for (clusterLabel in unassignedClusterLabels){
    X = clusters[Label == clusterLabel, X]
    Y = clusters[Label == clusterLabel, Y]
    clusterCentroid = colMeans(cbind(X, Y))
    unassignedClusterCentroids[Label == clusterLabel, X := clusterCentroid[1]]
    unassignedClusterCentroids[Label == clusterLabel, Y := clusterCentroid[2]]
    print("One iteration of for loop computing and storing unassigned cluster centroid.")
  }
  # Loop through assignedPoints, to see if any unassigned cluster centroids fall within Minor_Axis radius from assigned cluster centroid:
}




#trash:

    #control whether we are looping through prospective clusters:
    prospectiveClusterRemaining = TRUE
    # set predictedTrueClusterCentroid equal to the centroid of the cluster to which the point has been assigned:
    predictedTrueClusterCentroid = t(as.data.frame(colMeans(clusters[Label == point$cluster_ID, c(X, Y),])))
    # search within radius determined by point metadata to determine if an additional cluster is represented by the point:

    # if so, assign cluster to point:

    # "merge" clusters and compute new centroid of merged clusters:

    # ^ end inner loop: again search within radius determined by point metadata to determine if an additional cluster is represented by the point:






