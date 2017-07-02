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
source("/Users/seanhendryx/githublocal/assignPointsToClusters/assignPointsToClusters.R")


# Run test:

setwd("/Users/seanhendryx/DATA/Lidar/SRER/maxLeafAreaOctober2015/OPTICS_Param_Tests/study-area")

# read in clustered point cloud:
clusters = as.data.table(read.csv("OPTICS_clustered_points_eps_8.3_min_samples_150.csv"))
colnames(clusters)[1] = 'X'

# read in points:
points = as.data.table(read.csv("/Users/seanhendryx/DATA/SRERInSituData/SRER_Mesq_Tower_In_Situ_Allometry/inSituCoordinatesAndMeasurements.csv"))

points[,cluster_ID := NULL]
#first remove unnecessary points points:
validIDs = c(1:170)
validIDs = as.character(validIDs)
points = points[Sample_ID %in% validIDs,]

# RUN THESIS ALGORITHMS:
assignedPoints = assignPointsToClusters(points, clusters)

# Now run checkIfPointRepresentsMoreThanOneCluster
startTime = Sys.time()
clusters = checkIfPointRepresentsMoreThanOneCluster(assignedPoints, clusters)

write.csv(clusters, "clusters_assigned_to_points.csv")
endTime = Sys.time()
timeTaken = endTime - startTime
timeTaken

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






