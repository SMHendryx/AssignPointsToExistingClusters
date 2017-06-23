# Algorithms for assigning points in one dataset to clusters in another dataset 
# flagging unlikely correspondences based on a distance threshold.
# Authored by Sean Hendryx while working at the University of Arizona

#load packages:
library(data.table)

#define Euclidean distance function for working with two points/vectors in n-dim space:
eucDist <- function(x1, x2){
  sqrt(sum((x1 - x2) ^ 2))
} 

#algos:
####################################################################################################################################################################################
# vectorized:
assignPointsToClusters <- function(clusters, points, x_col_name = 'X', y_col_name = 'Y', cluster_id_col_name = 'Label'){
  # Algorithm assigns points, in a dataset $\bf{P}$, to the closet cluster in another dataset, $\bf{C}$,
  # flags unlikely correspondences based on distance threshold
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
  #points[,x_closestCentroid := double()]
  #points[,y_closestCentroid  := double()]
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

    print(paste0("closestMember: ", closestMember$Label))
    points[i,cluster_ID := closestMember$Label]
    print(paste0("closest cluster member distance = ", closestMember$distance_to_point))
    points[i,distance_to_closest_cluster_member := closestMember$distance_to_point]
    points[i, X_closest_cluster_member := closestMember$X]
    points[i, Y_closest_cluster_member := closestMember$Y]
  }
  return(points)
}





thresholdPoints <- function(points, thresholdType = "dominateMode"){
  ######################################################################################################################################################
  # Function computes threshold beyond which it is unlikely that the point corresponds to the cluster.  Expects a point set (data table) that has been assigned to clusters.
  ######################################################################################################################################################
  #first copy points and clusters to be modified locally (not by reference)
  points = copy(points)
  # Add boolean column indicating if the closest cluster is beyond threshold
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
  plot(density(closestPoints$distance_to_closest_cluster_member))
  
  if(thresholdType == "dominateMode"){
    threshold = dominateModeDist * buffer
  }else{
    threshold = meanDist * buffer
  }
  points[, closest_cluster_outside_threshold := (distance_to_closest_cluster_member > threshold)]
  #
  return(points)
}



}