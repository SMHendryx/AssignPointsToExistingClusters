# Matching_Points_And_Clusters
# Authored by Sean Hendryx while working at the University of Arizona

#load packages:
library(data.table)


#### Function definitions: ###############################################################################################################################################################################################################################################################

#define Euclidean distance function for working with two points/vectors in n-dim space:
euc.dist <- function(x1, x2){
  sqrt(sum((x1 - x2) ^ 2))
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
x <- runif(1000,0,100)
plot(density(x))
  abline(v=dmode(x))


#### End function definitions ###############################################################################################################################################################################################################################################################

# Run test:

setwd("/Users/seanhendryx/DATA/Lidar/SRER/maxLeafAreaOctober2015/OPTICS_Param_Tests/study-area")

# read in clustered point cloud:
clusters = as.data.table(read.csv("OPTICS_clustered_points_eps_8.3_min_samples_150.csv"))
colnames(clusters)[1] = 'X'

# read in points:
points = as.data.table(read.csv("/Users/seanhendryx/DATA/SRERInSituData/SRER_Mesq_Tower_In_Situ_Allometry/inSituCoordinatesAndMeasurements.csv"))

points[,cluster_ID := NULL]

assignedPoints = assignPointsToClusters(points, clusters)


################################################################################################################################################################################################################################################
#plot distance translation distribution:
library(ggplot2)
p = ggplot(assignedPoints, aes(x = distance_to_centroid)) + geom_density(fill = "#3ec09a", alpha = 0.5) + theme_bw() + labs(x = "Distance from Point to Cluster Centroid (m)", y = "Density")
p

# Plot assigned Points:


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

plotDT = clusters[Label != -1,]
plotDT = droplevels(plotDT)
ggp = ggplot(plotDT, aes(x = X, y = Y, color = Label)) + geom_point() + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 
ggp

###########################################################################################################################
#Making plot showing assingment of points to clusters:

#organize data to be rbinded:
#first remove unnecessary points from assignedPoints:
validIDs = c(1:170)
validIDs = as.character(validIDs)

assignedPoints = assignedPoints[Sample_ID %in% validIDs,]

ggp = ggplot() + geom_point(mapping = aes(x = X, y = Y, color = Label), data = plotDT) + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 

ggp = ggp + geom_point(mapping = aes(x = X, y = Y),data = assignedPoints, shape = 8)

# removing in situ points outside of study area:
maxX = max(plotDT[,X])
minX = min(plotDT[,X])
maxY = max(plotDT[,Y])
minY = min(plotDT[,Y])
assignedPoints = assignedPoints[X < maxX & X > minX & Y < maxY & Y > minY]







#trash:
files = list.files(pattern = "*.csv")

dir.create("Graphs")

for(file in files){
  clustered = as.data.table(read.csv(file))
  #labels = clustered[,Label]
  names(clustered)[1] = 'X'
  clustered[,Label := factor(Label)]
  ggp = ggplot(clustered, aes(x = X, y = Y, color = Label))
  ggp = ggp + geom_point() + theme_bw()
  name = substr(file, 1, nchar(file)-4)
  ggp
  ggsave(paste0("Graphs/",name, ".png"), device = 'png')
}























