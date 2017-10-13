# Algorithms for assigning points in one dataset to clusters in another dataset, flagging unlikely correspondences based on a distance threshold, 
# and "merging" clusters that have been over-segmented by the clustering algorithm.
# Ideally, if we have a dataset of points, where each point represents a cluster in another data set, there would be a one-to-one 
# correspondence between the two data sets, such that for each point there is one and only one cluster to which it corresponds and
# such that it is unmistakable which point corresponds to which cluster.  
# Though, this is not usually the case.  These algorithms present one approach to 
# 
# Authored by Sean Hendryx while working at the University of Arizona

#time it:
startTime = Sys.time()

#load packages:
library(data.table)
library(ggplot2)
source("/Users/seanhendryx/githublocal/assignPointsToClusters/APTEC.R")


#helper functions first:
printer <- function(string, variable){
  print(paste0(string, variable))
}

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
  
  # for checkIfPointRepresentsMoreThanOneCluster
  points[,x_closestCentroid := double()]
  points[,y_closestCentroid  := double()]
  
  # remove outliers coded as -1:
  clusters = clusters[Label != -1,]
  clusterLabels = unique(clusters[,cluster_id_col_name, with = FALSE])
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
  
  # Compute remaining unassigned cluster labels:
  unassignedClusterLabels = unique(clusters[is.na(assigned_to_point), Label])
  
  # Compute unassignedClusterCentroids:
  unassignedClusterCentroids = computeUnassignedClusterCentroids(clusters)

  for (unassignedClusterLabel in unassignedClusterLabels){
    # test if point falls within minor axis circle from assigned cluster centroid:
    #x, center_x, y, center_y, radius
    print(paste0("Testing unassigned cluster: ", unassignedClusterLabel))
    x = unassignedClusterCentroids[Label == unassignedClusterLabel, X]
    y = unassignedClusterCentroids[Label == unassignedClusterLabel, Y]

    #printer("x: ", x)
    #printer("center_x: ", center_x)
    #printer("y: ", y)
    #printer("center_y: ", center_y)
    #printer("radius: ", radius)
    if (testIfPointWithinCircle(x = x, center_x = center_x, y = y, center_y = center_y, radius = radius)){
      print(paste0("Cluster centroid falls within radius."))
      #if unassigned cluster centroid within minor_axis radius of assigned centroid,
      # assign point to cluster:
      # ERROR HERE: i am here:
      clusters[Label == unassignedClusterLabel, assigned_to_point := assignedPoints[Sample_ID == pointID, Sample_ID]]
      
      # compute newPredictedCentroid from clusters:
      X = clusters[assigned_to_point == pointID, X]
      Y = clusters[assigned_to_point == pointID, Y]
      newPredictedCentroid = colMeans(cbind(X, Y))
      print(paste0("New predicted centroid: ", newPredictedCentroid))

      #Recursive call:
      print("Starting recursive call to testAndMergeClustersRecursively:")
      testAndMergeClustersRecursively(newPredictedCentroid, pointID, assignedPoints, clusters)
    }# end if (testIfPointWithinCircle)
  }# end for unassignedClusterLabel in unassignedClusterLabels
}

checkIfPointRepresentsMoreThanOneCluster <- function(assignedPoints, clusters){
  # Function resolves over-segmentation of clusters by determining if any other clusters should be assigned to the points in assignedPoints
  # based on information held in the point (metadata) and, if so,
  # assigns the cluster(s) to the point.
  # Assumes that metadata is at the same scale of clusters
  # Updates the clusters data.table object by reference.
  
  # First, remove clusters outliers coded as -1:
  clusters = clusters[Label != -1,]
  # since factor, Label -1 still exists as level, so:
  clusters[,Label := droplevels(Label)]

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
  rm(uniquelyAssignedPointIDs)
  rm(uniquelyAssignedPointIDs_i)
  uniquelyAssignedPointIDs = vector(mode = "character")
  assignedClusterIDs = unique(assignedPoints$cluster_ID)
  i = 1
  for (id in assignedClusterIDs){
    if(nrow(assignedPoints[cluster_ID == id]) == 1) {
      uniquelyAssignedPointIDs_i = assignedPoints[cluster_ID == id, Sample_ID]
      #DT[["region"]]
      #i = i + 1 
      uniquelyAssignedPointIDs = c(uniquelyAssignedPointIDs, uniquelyAssignedPointIDs_i)
      #printer("uniquelyAssignedPointIDs: ", uniquelyAssignedPointIDs)
      print(uniquelyAssignedPointIDs)
    } 
  }
  
  #I am here
  for(pointID in uniquelyAssignedPointIDs){
    print(paste0("Testing if any other clusters belong to point: ", pointID))
    # compute predictedCentroid of true cluster from clusters (i.e., compute the centroid of the points in cluster datatable that are assigned to point with pointID):
    X = clusters[assigned_to_point == pointID, X]
    Y = clusters[assigned_to_point == pointID, Y]
    predictedCentroid = colMeans(cbind(X, Y))

    testAndMergeClustersRecursively(predictedCentroid = predictedCentroid, pointID = pointID, assignedPoints = assignedPoints, clusters = clusters)
  }
}


#### End function definitions ###############################################################################################################################################################################################################################################################

# Run test:

setwd("/Users/seanhendryx/DATA/Lidar/SRER/maxLeafAreaOctober2015/OPTICS_Param_Tests/study-area")

# read in clustered point cloud:
clusters = as.data.table(read.csv("OPTICS_clustered_points_eps_8.3_min_samples_150.csv"))
colnames(clusters)[1] = 'X'

# read in points:
points = as.data.table(read.csv("/Users/seanhendryx/DATA/SRERInSituData/SRER_Mesq_Tower_In_Situ_Allometry/inSituCoordinatesAndMeasurements.csv"))

points[,cluster_ID := NULL]

# RUN THESIS ALGORITHMS:
assignedPoints = assignPointsToClusters(points, clusters)

# Now run checkIfPointRepresentsMoreThanOneCluster
startTime = Sys.time()
checkIfPointRepresentsMoreThanOneCluster(assignedPoints, clusters)
endTime = Sys.time()
timeTaken = endTime - startTime
print(timeTaken)

################################################################################################################################################################################################################################################
#####  PLOTS ##################################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################

# Plot assigned & threshed Points over clusters:
clusters$Label = factor(clusters$Label)
# make qualitative color palette:
# 82 "color blind friendly" colors from: http://tools.medialab.sciences-po.fr/iwanthue/
# with outliers set to black
cbf = c("#000000", "#be408c",
  "#4cdc8b",
  "#b1379e",
  "#90d15e",
  "#442986",
  "#dab528",
  "#577ceb",
  "#cba815",
  "#424cad",
  "#acbd3d",
  "#745bc4",
  "#7bcf6e",
  "#863c9f",
  "#4eaa4c",
  "#e768c5",
  "#669b2c",
  "#9e7ee9",
  "#2e7c23",
  "#c180e2",
  "#a6bf55",
  "#6d1b66",
  "#37d8b0",
  "#c42a6d",
  "#5dba6f",
  "#8d3f90",
  "#af9a23",
  "#6d7ddb",
  "#e7af45",
  "#468ae0",
  "#dd8026",
  "#4e62aa",
  "#c1851b",
  "#3d3072",
  "#bdb553",
  "#835fb5",
  "#70851b",
  "#e592e6",
  "#255719",
  "#b765b8",
  "#3eac74",
  "#992963",
  "#77daa9",
  "#ab2540",
  "#36dee6",
  "#cd3e43",
  "#3aad8c",
  "#e25968",
  "#458541",
  "#db81c4",
  "#516f1d",
  "#c093db",
  "#817614",
  "#7199e0",
  "#a54909",
  "#894f8e",
  "#9fc069",
  "#6b1740",
  "#8cbf79",
  "#d95987",
  "#b9b567",
  "#97436d",
  "#e0b75e",
  "#de7bae",
  "#818035",
  "#d04c6c",
  "#b18b34",
  "#e67d9f",
  "#a06919",
  "#822131",
  "#d0a865",
  "#7d2716",
  "#e29249",
  "#c76674",
  "#80591b",
  "#e77162",
  "#9e3119",
  "#e1925d",
  "#d26c69",
  "#d16c2f",
  "#c46d53",
  "#e26a4a",
  "#aa612f")

# Plot only those points inside threshold:
#organize data to be rbinded:
#first remove unnecessary points from assigned and thresholded Points:
validIDs = c(1:170)
validIDs = as.character(validIDs)
assignedPoints = assignedPoints[Sample_ID %in% validIDs,]

#remove outliers (coded -1) in clusters data.table:
plotDT = clusters[Label != -1,]
plotDT = droplevels(plotDT)

# removing in situ points outside of study area:
maxX = max(plotDT[,X])
minX = min(plotDT[,X])
maxY = max(plotDT[,Y])
minY = min(plotDT[,Y])
assignedPoints = assignedPoints[X < maxX & X > minX & Y < maxY & Y > minY]


renderStartTime = Sys.time()
ggp = ggplot() + geom_point(mapping = aes(x = X, y = Y, color = factor(Label)), data = plotDT, size = .75) + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 
#testing adding assigned_to_point column to clusters:
#ggp = ggplot() + geom_point(mapping = aes(x = X, y = Y, color = factor(assigned_to_point)), data = plotDT, size = .75) + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 

ggp = ggp + geom_point(data = assignedPoints[closest_cluster_outside_threshold == FALSE,], mapping = aes(x = X, y = Y), shape = 8)
ggp = ggp + geom_point(data = assignedPoints[closest_cluster_outside_threshold == FALSE,], mapping = aes(x = X_closest_cluster_centroid, y = Y_closest_cluster_centroid), shape = 13)
ggp

endTime = Sys.time()
renderTimeTaken = endTime - renderStartTime
timeTaken = endTime - startTime
print("Time taken to render graph: ")
print(renderTimeTaken)

print("Total Time Taken: ")
print(timeTaken)














######################################################################################################################################################################################################################################################

#plot distance translation distribution:
p = ggplot(testAssignedPoints, aes(x = distance_to_closest_cluster_member)) + geom_density(fill = "#3ec09a", alpha = 0.5) + theme_bw() + labs(x = "Distance from Point to Closest Cluster Member (m)", y = "Density")
p

######################################################################################################################################################################################################################################################

#plot all clusters:
ggp = ggplot(plotDT, aes(x = X, y = Y, color = Label)) + geom_point() + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 
ggp

######################################################################################################################################################################################################################################################
#Making plot showing assignment of points to clusters:

#organize data to be rbinded:
#first remove unnecessary points from assignedPoints:
validIDs = c(1:170)
validIDs = as.character(validIDs)

assignedPoints = assignedPoints[Sample_ID %in% validIDs,]

ggp = ggplot() + geom_point(mapping = aes(x = X, y = Y, color = Label), data = plotDT) + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 

ggp = ggp + geom_point(mapping = aes(x = X, y = Y),data = assignedPoints, shape = 8)

# removing in situ points outside of study area:
# NOW COMPLETED BEFORE RUNNING assignPointsToClusters()
#maxX = max(plotDT[,X])
#minX = min(plotDT[,X])
#maxY = max(plotDT[,Y])
#minY = min(plotDT[,Y])
#assignedPoints = assignedPoints[X < maxX & X > minX & Y < maxY & Y > minY]




###########################################################################################################################


p = ggplot(threshedPoints[closest_cluster_outside_threshold == FALSE,], aes(x = distance_to_centroid)) + geom_density(fill = "gray41", alpha = 0.5) + theme_bw() + labs(x = "Distance from Point to Cluster Centroid (m)", y = "Density")
p






