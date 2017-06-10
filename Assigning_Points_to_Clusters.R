# Matching_Points_And_Clusters
# Authored by Sean Hendryx while working at the University of Arizona

#load packages:
library(data.table)

#define Euclidean distance function for working with two points/vectors in n-dim space:
euc.dist <- function(x1, x2){
  sqrt(sum((x1 - x2) ^ 2))
} 
assignPointsToClustersWithThreshold <- function(points, clusters, x_col_name = 'X', y_col_name = 'Y', cluster_id_col_name = 'Label', compute_centroids = TRUE){
  # this algorithm assigns point values to clusters.  For example, 
  # if we have a matrix of point coordinates, each of which should represent
  # a cluster in some way.
  # if outliers are coded as -1 in cluster_id_col_name, they will be assumed to not be clusters
  # :Param points: data.table object with columns 'X' and 'Y'
  # :Param clusters: data.table object with columns 'X', 'Y', and 'Label'
  #check if column already exists:
  if(any(names(points) == "cluster_ID")){
    stop("cluster_ID already exists in points.  points should not include cluster ids prior to running assignPointsToClusters function.")
  }
  #if doesn't exist, add:
  points[,cluster_ID := integer()]
  points[,x_closestCentroid := double()]
  points[,y_closestCentroid  := double()]
  clusterLabels = unique(clusters[,cluster_id_col_name, with = FALSE])
  # remove outliers coded as -1:
  clusterLabels = clusterLabels[eval(parse(text=cluster_id_col_name)) != -1,]
  #print(paste0("clusterLabels", clusterLabels))
  for(i in seq(nrow(points))){
    print(paste0("Looping through points, on point: ", i, "\n"))
    position = points[i, c(x_col_name, y_col_name), with = FALSE]
    minDist = Inf
    #for each point, find the closest cluster centroid:
    for(cluster in clusterLabels[,eval(parse(text=cluster_id_col_name))]){
      #print(paste0("Looping through clusters, on cluster: ", cluster, "\n"))
      if(compute_centroids){
        #print("cluster: \n")
        #print(cluster)
        #print("\n")
        #print(clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE])
        cluster_centroid = t(as.data.frame(colMeans(clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE])))
        #print(paste0("cluster_centroid ", cluster_centroid, "\n"))
      }else{
        cluster_centroid = clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE]
      }
      distance = euc.dist(position, cluster_centroid)
      #print(paste0("distance from point to cluster centroid = ", distance, "\n"))
      if(distance < minDist){
        closestCluster = cluster
        closestCentroid_x = cluster_centroid[1]
        closestCentroid_y = cluster_centroid[2]
        minDist = distance
      }
    }
    #assign closest cluster to point:
    # because we only loop through each point once, but loop through clusters for each point, assigning in this fashion effectively
    # assigns the points to the clusters because one cluster can be assigned multiple points while each point can be assigned only one cluster
    print(paste0("closest cluster = ", closestCluster))
    print(paste0("closest cluster distance (minDist) = ", minDist))
    points[i,cluster_ID := closestCluster]
    points[i,x_closestCentroid := closestCentroid_x]
    points[i,y_closestCentroid := closestCentroid_y]
    points[i,distance_to_centroid := minDist]
  }
  ######################################################################################################################################################
  # Now compute threshold beyond which it is unlikely that the point corresponds to the cluster:
  #find the distance of the closest in situ coordinate (i.e. "point" in points) to each cluster centroid:
  cluster_IDs = unique(points[,cluster_ID])
  #instantiate list of distances from each cluster to it’s closest point (which we will take the mean of):
  cluster_distances = c()
  for (cluster in cluster_IDs){
    #create index i for indexing into cluster_distances list:
    i = 1
    #actually, just use min(points_assigned_to_cluster[,c(distance_to_centroid, Sample_ID)])
    print(min(points[cluster_ID == clusterc(distance_to_centroid, Sample_ID)]))
    #here
    } 
  }

  return(points)
}

assignPointsToClusters <- function(points, clusters, x_col_name = 'X', y_col_name = 'Y', cluster_id_col_name = 'Label', compute_centroids = TRUE){
  # this algorithm assigns point values to clusters.  For example, 
  # if we have a matrix of point coordinates, each of which should represent
  # a cluster in some way.
  # if outliers are coded as -1 in cluster_id_col_name, they will be assumed to not be clusters
  # :Param points: data.table object with columns 'X' and 'Y'
  # :Param clusters: data.table object with columns 'X', 'Y', and 'Label'
  #check if column already exists:
  if(any(names(points) == "cluster_ID")){
    stop("cluster_ID already exists in points.  points should not include cluster ids prior to running assignPointsToClusters function.")
  }
  #if doesn't exist, add:
  points[,cluster_ID := integer()]
  points[,x_closestCentroid := double()]
  points[,y_closestCentroid  := double()]
  clusterLabels = unique(clusters[,cluster_id_col_name, with = FALSE])
  # remove outliers coded as -1:
  clusterLabels = clusterLabels[eval(parse(text=cluster_id_col_name)) != -1,]
  #print(paste0("clusterLabels", clusterLabels))
  for(i in seq(nrow(points))){
    print(paste0("Looping through points, on point: ", i, "\n"))
    position = points[i, c(x_col_name, y_col_name), with = FALSE]
    minDist = Inf
    #for each point, find the closest cluster centroid:
    for(cluster in clusterLabels[,eval(parse(text=cluster_id_col_name))]){
      #print(paste0("Looping through clusters, on cluster: ", cluster, "\n"))
      if(compute_centroids){
        #print("cluster: \n")
        #print(cluster)
        #print("\n")
        #print(clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE])
        cluster_centroid = t(as.data.frame(colMeans(clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE])))
        #print(paste0("cluster_centroid ", cluster_centroid, "\n"))
      }else{
        cluster_centroid = clusters[eval(parse(text=cluster_id_col_name)) == cluster, c(x_col_name, y_col_name), with = FALSE]
      }
      distance = euc.dist(position, cluster_centroid)
      #print(paste0("distance from point to cluster centroid = ", distance, "\n"))
      if(distance < minDist){
        closestCluster = cluster
        closestCentroid_x = cluster_centroid[1]
        closestCentroid_y = cluster_centroid[2]
        minDist = distance
      }
    }
    #assign closest cluster to point:
    # because we only loop through each point once, but loop through clusters for each point, assigning in this fashion effectively
    # assigns the points to the clusters because one cluster can be assigned multiple points while each point can be assigned only one cluster
    print(paste0("closest cluster = ", closestCluster))
    print(paste0("closest cluster distance (minDist) = ", minDist))
    points[i,cluster_ID := closestCluster]
    points[i,x_closestCentroid := closestCentroid_x]
    points[i,y_closestCentroid := closestCentroid_y]
    points[i,distance_to_centroid := minDist]
  }
  return(points)
}





















