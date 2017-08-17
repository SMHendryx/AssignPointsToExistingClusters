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


#generate clusters with true cluster identity (cluster) and a made-up clustering (clusterHat):
c1Centroid = c(2,2)
c2Centroid = c(4,2)
c3Centroid = c(6,3)

