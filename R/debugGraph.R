#load packages:
library(data.table)
library(ggplot2)

setwd("/Users/seanhendryx/DATA/Lidar/SRER/maxLeafAreaOctober2015/OPTICS_Param_Tests/study-area")

# read in clustered point cloud:
clusters = as.data.table(read.csv("OPTICS_clustered_points_eps_8.3_min_samples_150.csv"))
colnames(clusters)[1] = 'X'

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

ggp = ggplot() + geom_point(mapping = aes(x = X, y = Y, color = factor(Label)), data = plotDT, size = .75) + theme_bw() + theme(legend.position="none") + scale_colour_manual(values = cbf) 


# debug scraps:
# for warning: if (dx > radius) { ... :
#		  the condition has length > 1 and only the first element will be used
x= 512617.073694242
center_x= 512628.276596064
y= 3520666.88610194
center_y= 3520675.59270331
radius= 4.6

dx = abs(x-center_x)
dy = abs(y-center_y)
if(dx > radius){
print(FALSE)
}
if (dy > radius){
print(FALSE)
}
if (dx + dy <= radius){
print(TRUE)
}
if (dx^2 + dy^2 <= radius^2){
print(TRUE)
} else {
print(FALSE)
}
