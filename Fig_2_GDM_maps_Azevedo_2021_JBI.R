# Scripts for mapping observed and predicted turnover from Azevedo et al., 2021
# Here you you'll find how to produce Figure 2
# Author: Josué A. R. Azevedo
# Department of biological and environmental sciences - Uni of Gothenburg
# GGBC Sweden
# INPA AM Brazil
# email: josueanderson21@gmail.com

# Loading libraries
library(viewport)
library(grid)
library(gridExtra)
library(RStoolbox)
library(sf)
library(ggplot2)
library(ggpubr)
library(cluster)
library(rgeos)
library(dendextend)

# Set your working directory to the folder downloaded from GitHub
setwd("")

# Read in cerrado_layers.rds downloaded from GitHub ####

cerrado_layers = readRDS("cerrado_layers.rds")

# Within this list, you will find the following objects in this order (all spatial objects are on an equal-area Behrmann projection):

grid_cerr_three_or_more_an = cerrado_layers[[1]][[1]]  # grid with localities for amphibians

grid_cerr_three_or_more_sq = cerrado_layers[[1]][[2]]  # grid with localities for reptiles

grid_todos = cerrado_layers[[1]][[3]]  # grid with localities for both groups

elevation_c = cerrado_layers[[1]][[4]]  # raster with elevation values

slope_roughness_c = cerrado_layers[[1]][[4]] # raster with slope values

relief_roughness_c = cerrado_layers[[1]][[6]] # raster with relief roughness values

precip_rast = cerrado_layers[[1]][[7]]  # raster with precipitation variables (CHELSA)

temp_rast = cerrado_layers[[1]][[8]]  # raster with temperature variables (CHELSA)

soil_raster = cerrado_layers[[1]][[9]]  # raster with soil variables (%sand, %clay, %coarse fragments; soilgrids.org)

All_sq_phylo.beta.sim_FULLY = cerrado_layers[[2]]  # matrix of phylogenetic turnover for reptiles

Endemic_sq_phylo.beta.sim_FULLY = cerrado_layers[[3]] # matrix of phylogenetic turnover for endemic reptiles 

Endemic_an_phylo.beta.sim_Full = cerrado_layers[[4]]  # matrix of phylogenetic turnover for endemic amphibians

All_an_phylo.beta.sim_Full = cerrado_layers[[5]]  # matrix of phylogenetic turnover for amphibians

all_an.betapair = cerrado_layers[[6]]  # matrix of taxonomic turnover for amphibians

END_an.beta.sim = cerrado_layers[[7]]  # matrix of taxonomic turnover for endemic amphibians

all_sq.beta.sim = cerrado_layers[[8]] # matrix of taxonomic turnover for reptiles

END_sq.beta.sim = cerrado_layers[[9]]  # matrix of taxonomic turnover for endemic amphibians

pd_sd_sq = cerrado_layers[[10]] # standard deviation of PD per location across a sample of 100 phylogenies for reptiles

pd_sd_an = cerrado_layers[[11]] # standard deviation of PD per location across a sample of 100 phylogenies for amphibians

cerrado = cerrado_layers[[12]] # Shapefile delimiting the Cerrado Savanna borders

worldMap = cerrado_layers[[13]] # Multipolygon sf object delimiting the continental borders without small islands

######

# Defining Behrmann equal area projection
behr <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'
# Transforming spatial polygon dataframe into a sf object
worldMap= st_as_sf(worldMap)
map_data_sf = st_transform(worldMap, "+proj=cea +lat_ts=30 +ellps=WGS84 ")
map_data_sf_s <- rmapshaper::ms_simplify(input = as(map_data_sf, 'Spatial')) %>% st_as_sf()

###########################################################
######## Figure 2 in the paper ############################
###########################################################

## Plotting observed clustering of phylogentic turnover (locs) ####
#
# General settings####
k = 12 # Define number of clusters for UPGMA clustering
size_pts= 3 # Size of the point localities for plotting
cerrado_buf = raster::buffer(cerrado,500000) # Buffer around Cerrado layer for delimiting the extent
# Refining limits for plotting (particular to my dataset)
min_x=extent(cerrado_buf)[1]-000000
max_x=extent(cerrado_buf)[2]-000000
min_y=extent(cerrado_buf)[3]+600000
max_y=extent(cerrado_buf)[4]+300000
cerrado_sf = st_as_sf(cerrado) # Transform from SpatialPolygon to sf object
######

#
# Fig. 1a - Turnover for reptiles ####
expression_title = paste0("(a) Reptiles - ", "Observed ", "Phylo-turnover")
# PhyloTurnover for reptiles
turnover_matrix =  All_sq_phylo.beta.sim_FULLY
# Getting the UPGMA clustering ####
tree = agnes(turnover_matrix, method = "average")
cluster_membership<- cutree(tree, k = k)

# Transforming grid centroids with samples in points for plotting 
grid_todos_end= grid_todos[grid_todos@data$layer %in% colnames(as.matrix(turnover_matrix)),]
pts_tods = gCentroid(grid_todos_end, byid=TRUE)

# Assigning cluster membership for each locality
pts_df = as.data.frame(cbind(pts_tods@coords,cluster_membership))
head(pts_df)
pts_df$cluster_membership = as.factor(cluster_membership)
#####
## Preparing a UPGMA dendrogram for plotting ####
dend <- as.dendrogram(tree)

# Checking the order of the colnames ###
colnames(as.matrix(turnover_matrix))
names(cluster_membership)
labels(dend)
names(cluster_membership) == labels(dend) # Not in the same order

# Change cluster membership label order according to the dend label order
cluster_membership2=cluster_membership[order(match(names(cluster_membership), labels(dend)))]
head(cluster_membership2)

# Check the order again
names(cluster_membership2) == labels(dend)

dend2=branches_attr_by_clusters(dend,cluster_membership2)

# This is to visualize to which branches the clusters are assigned
dend4= dend2
labels(dend4) = cluster_membership2
dend4 %>% set("branches_lwd",2.5) %>% ladderize(TRUE)  %>%  plot(horiz  = FALSE,labels=FALSE,axes=FALSE)

# to get the right order from of clusters for assigning colors
unique(cluster_membership2)

# Order our point localities (grid cells centroids) according to the dend order
pts_df= pts_df[order(factor(pts_df$cluster_membership,levels=c(unique(cluster_membership2)))),]

# Assigning colors for the dendrogram branches. This demands a bit of manual work to assing divergent colors to increasingly distant clusters if you are using your own dataset.
dend3 <- color_branches(dend2,clusters = cluster_membership2,  col=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Create a column in pts_df to define colors for plotting. Same order as the dendrogram.
pts_df$col= plyr::mapvalues(pts_df$cluster_membership, from=unique(cluster_membership2), to=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Checking the order of the colors in the objects
unique(cluster_membership2)
unique(get_leaves_branches_col(dend3))

# Produce a ggplot object to be added the final map
par(bg="transparent")
dude= dend3 %>% set("branches_lwd",.4) %>% set("labels", NA)  %>% ggplot(horiz  = TRUE,labels=FALSE,axes=FALSE)
ggplot(as.ggdend(dend3, horiz  = TRUE,labels=FALSE,axes=FALSE))
#####
# Plotting itself #####
expression_title = paste0("(a) Reptiles- Observed Phylogenetic Turnover")
reptiles_plot= ggplot(map_data_sf)+
  geom_sf(data=map_data_sf,fill="grey60",
          col="grey30",
          size=0.1)+
  geom_sf(data=cerrado_sf,
          colour="grey",
          fill="white",
          size=.25)+
  geom_point(data=pts_df,aes(x=x, y=y),shape=21, fill = pts_df$col,size=size_pts)+
  geom_sf(data=map_data_sf,
          colour="grey30",
          fill=NA,size=0.1)+
  coord_sf(xlim=c(min_x,max_x),ylim = c(min_y,max_y))+
  scale_x_continuous(breaks = seq(-70, -40, by = 10))+
  scale_y_continuous(breaks = seq(-30, 10, by = 10))+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'white',size = .5),
        panel.background = element_rect(fill = 'grey95'))+
  theme(panel.border = element_rect(colour = NA, fill=NA, size=.5))+
  theme(axis.text=element_text(size=8))+
  theme(axis.line = element_line(colour = 'grey20', size = .1))+
  theme(legend.background=element_rect(fill = NA, colour = NA),
        legend.title = element_text(color="grey20",family="Times"),
        legend.text  = element_text(color="grey20",family="Times"),
        legend.key.height=unit(0.04,"npc"),
        legend.key.width=unit(0.01,"npc"),
        legend.position = c(1.0001,.75),
        legend.justification = 'left')+
  geom_rect(xmin = -Inf, xmax = -5050000,   ymin = 370000,    ymax = -1400000, fill = "grey60") +
  annotate(geom="text", x=-5800000, y= 295000, label="UPGMA clustering",
           color="white",size = 3.5,family="Times")+
  ggtitle(expression_title)+ 
  theme(plot.title = element_text(size=10,face="plain",hjust = 0,family="Times"))+
  theme(plot.title = element_text(face="plain"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
#reptiles_plot # The dendrogram is added only at the final plot due to adjustiments in the positions of each figure.

#####

#####
# Fig. 1B - Observed phylogenetic Turnover AMPHIBIANS ####
# Steps are better described for Fig 1A
turnover_matrix = All_an_phylo.beta.sim_Full
k = 12
expression_title = paste0("(b) Amphibians - ", "Observed ", "Phylo-turnover")
#####
# Getting the cluster ####
tree = agnes(turnover_matrix, method = "average")
cluster_membership<- cutree(tree, k = k)

# Transforming in points for plotting 
grid_todos_end= grid_todos[grid_todos@data$layer %in% colnames(as.matrix(turnover_matrix)),]
pts_tods = gCentroid(grid_todos_end, byid=TRUE)
pts_df = as.data.frame(cbind(pts_tods@coords,cluster_membership))
head(pts_df)
pts_df$cluster_membership = as.factor(cluster_membership)
#####
## Getting a dendrogram ####
dend <- as.dendrogram(tree)

# Checking the order of the colnames ###
head(colnames(as.matrix(turnover_matrix)))
head(cluster_membership)
order.dendrogram(dend)

# Changing cluster membership according to the dend order
cluster_membership2=cluster_membership[order(match(names(cluster_membership), labels(dend)))]

dend2=branches_attr_by_clusters(dend,cluster_membership2)

# This is to visualize to which branches are the clusters assigned
dend4= dend2
labels(dend4) = cluster_membership2
dend4 %>% set("branches_lwd",2.5) %>% ladderize(TRUE)  %>%  plot(horiz  = FALSE,labels=FALSE,axes=FALSE)

# to get the right order from the phylogenetic turnover mapping
unique(cluster_membership2)
unique(get_leaves_branches_col(dend3))

# Assigning colors to the localities.
pts_df$col= plyr::mapvalues(pts_df$cluster_membership, from=unique(cluster_membership2), to=c(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue")))

# Color the dendrogram
dend3 <- color_branches(dend2,clusters = cluster_membership2,  col=c(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue")))

# Testning the order of the colors
unique(cluster_membership2)
unique(get_leaves_branches_col(dend3))

# Produce a ggplot for the dendrogram
par(bg="transparent")
dude_an= dend3 %>% set("branches_lwd",.4) %>% set("labels", NA)  %>% ggplot(horiz  = TRUE,labels=FALSE,axes=FALSE)
ggdend(dend3, horiz  = TRUE,labels=FALSE,axes=FALSE)
#####
# Plotting itself #####
expression_title = paste0("(b) Amphibians- Observed Phylogenetic Turnover")
anf_plot= ggplot(map_data_sf)+
  geom_sf(data=map_data_sf,fill="grey60",
          col="grey30",
          size=0.1)+
  geom_sf(data=cerrado_sf,
          colour="grey",
          fill="white",
          size=.25)+
  geom_point(aes(x=x, y=y),shape=21, fill = pts_df$col,data=pts_df,size=size_pts)+
  geom_sf(data=map_data_sf,
          colour="grey30",
          fill=NA,size=0.1)+
  coord_sf(xlim=c(min_x,max_x),ylim = c(min_y,max_y))+
  scale_x_continuous(breaks = seq(-70, -40, by = 10))+
  scale_y_continuous(breaks = seq(-30, 10, by = 10))+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'white',size = .5),
        panel.background = element_rect(fill = 'grey95'))+
  theme(panel.border = element_rect(colour = NA, fill=NA, size=.5))+
  theme(axis.text=element_text(size=8))+
  theme(axis.line = element_line(colour = 'grey20', size = .1))+
  theme(legend.background=element_rect(fill = NA, colour = NA),
        legend.title = element_text(color="grey20",family="Times"),
        legend.text  = element_text(color="grey20",family="Times"),
        legend.key.height=unit(0.04,"npc"),
        legend.key.width=unit(0.01,"npc"),
        legend.position = c(1.0001,.75),
        legend.justification = 'left')+
  ggtitle(expression_title)+ 
  geom_rect(xmin = -Inf, xmax = -5050000,   ymin = 370000,    ymax = -1400000, fill = "grey60") +
  annotate(geom="text", x=-5800000, y= 295000, label="UPGMA clustering",
           color="white",size = 3.5,family="Times")+
  theme(plot.title = element_text(size=10,face="plain",hjust = 0,family="Times"))+
  theme(plot.title = element_text(face="plain"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
#anf_plot####
#####

#####
# Fig 1C) PREDICETED Clustering of phylo turnover for reptiles ####
prep= c("bio5","bio15") # Variables selected in the best model for reptiles
# I will run gdm for the best model. The turnover matrix for all reptiles is All_sq_phylo.beta.sim_FULLY. Remember to load the GDM_DF_function from the script phylo_turnover_azevedo_2021.R
end_rep_gdm_table = GDM_DF_function(turnover_matrix = All_sq_phylo.beta.sim_FULLY,
                                    weights_mat1 = pd_sd_sq,
                                    best_model = TRUE,
                                    var_to_final_model= prep)
gdm_end_rep <- gdm::gdm(end_rep_gdm_table, geo=TRUE)

# Getting rasters of all variables
preds=stack(elevation_c, slope_roughness_c,relief_roughness_c , precip_rast,temp_rast,soil_raster)

# Transforming rasters according to GMD predictions
rastTrans_end_rep <- gdm::gdm.transform(gdm_end_rep, subset(preds,prep))

# Masking the rasters to the Cerrado borders
rastTrans_end_rep<- raster::mask(rastTrans_end_rep,cerrado)
plot(rastTrans_end_rep)

# Getting point values from the transformed raster
sampling_pts = rasterToPoints(rastTrans_end_rep[[1]])
sampling_pts=sampling_pts[,-3]
where_about = raster::extract(rastTrans_end_rep,sampling_pts)
sampling_pts_2 = cbind(sampling_pts,where_about)
sampling_pts_3 = na.exclude(sampling_pts_2)
attr(sampling_pts_3, "ATT") <- NULL
head(sampling_pts_3)

# Preparing dataframe of gdm transformed variables for UPGMA
for_agnes = sampling_pts_3[,-c(1:2)]
tree_pred_rep = agnes(for_agnes, method = "average") # This one takes a while
cluster_membership2=cutree(tree_pred_rep, k = k)

# Getting cluster memberships delimited above and assign it to each point
pts_to_plot= as.data.frame(cbind(sampling_pts_3[,c(1:2)],cluster_membership2))
pts_to_plot$cluster_membership2 = cluster_membership2
coordinates(pts_to_plot)=pts_to_plot[,c(1:2)]
pred_rast = raster::rasterize(pts_to_plot,rastTrans_end_rep,field=pts_to_plot@data$cluster_membership)
pred_pol = st_as_sf(rasterToPolygons(raster::rasterize(pts_to_plot,rastTrans_end_rep,field=pts_to_plot@data$cluster_membership2),dissolve = TRUE))
pts_to_plot_df = pts_to_plot@data
plot(pred_rast)
plot(pred_pol, add=TRUE)
#####
## Getting a dendrogram (steps are better described in Fig. 1A) ####
dend <- as.dendrogram(tree_pred_rep)

# Checking the order of the colnames ###
head(colnames(as.matrix(All_sq_phylo.beta.sim_FULLY)))
head(cluster_membership)
order.dendrogram(dend)

# Changing cluster membership according to the dend order
cluster_membership3=cluster_membership2[order(match(names(cluster_membership2), labels(dend)))]

dend2=branches_attr_by_clusters(dend,cluster_membership3)

# This is to visualize to which branches are the clusters assigned
dend4= dend2
labels(dend4) = cluster_membership3
#dend4 %>% set("branches_lwd",2.5) %>% ladderize(TRUE)  %>%  plot(horiz  = FALSE,labels=FALSE,axes=FALSE)

# to get the right order from the phylogenetic turnover mapping
unique(cluster_membership3)
unique(get_leaves_branches_col(dend3))

# To map values in the the points (will become a raster in ggplot) 
pts_to_plot_df$col= plyr::mapvalues(pts_to_plot_df$cluster_membership2, from=unique(cluster_membership3), to=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Color the dendrogram with the same colors above
dend3 <- color_branches(dend2,clusters = cluster_membership3,  col=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Testning the order of the colors
unique(pts_to_plot_df$cluster_membership)

unique(cluster_membership2)

# For ggplot in the same order above
cols=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue")

# Produce a ggplot for the dendrogram
par(bg="transparent")
dude_pred= dend3 %>% set("branches_lwd",.4) %>% set("labels", NA)  %>% ggplot(horiz  = TRUE,labels=FALSE,axes=FALSE)
ggdend(dend3, horiz  = TRUE,labels=FALSE,axes=FALSE)

# Plotting itself #####
expression_title = paste0("(c) Reptiles- Predicted Phylogenetic Turnover")
pred_reptiles_plot= ggplot(map_data_sf)+
  geom_sf(data=map_data_sf,fill="grey60",
          col="grey30",
          size=0.1)+
  geom_raster(data =pts_to_plot_df,
              aes(x=x,
                  y=y),
              fill= pts_to_plot_df$col)+
  geom_sf(data=map_data_sf,
          colour="grey30",
          fill=NA,size=0.1)+
  geom_sf(data=pred_pol,colour="black",fill=NA,size=0.2)+
  coord_sf(xlim=c(min_x,max_x),ylim = c(min_y,max_y))+
  scale_x_continuous(breaks = seq(-70, -40, by = 10))+
  scale_y_continuous(breaks = seq(-30, 10, by = 10))+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'white',size = .5),
        panel.background = element_rect(fill = 'grey95'))+
  theme(panel.border = element_rect(colour = NA, fill=NA, size=.5))+
  theme(axis.text=element_text(size=8))+
  theme(axis.line = element_line(colour = 'grey20', size = .1))+
  theme(legend.background=element_rect(fill = NA, colour = NA),
        legend.title = element_text(color="grey20",family="Times"),
        legend.text  = element_text(color="grey20",family="Times"),
        legend.key.height=unit(0.04,"npc"),
        legend.key.width=unit(0.01,"npc"),
        legend.position = c(1.0001,.75),
        legend.justification = 'left')+
  ggtitle(expression_title)+ 
  geom_rect(xmin = -Inf, xmax = -5050000,   ymin = 370000,    ymax = -1400000, fill = "grey60") +
  annotate(geom="text", x=-5800000, y= 295000, label="UPGMA clustering",
           color="white",size = 3.5,family="Times")+
  theme(plot.title = element_text(size=10,face="plain",hjust = 0,family="Times"))+
  theme(plot.title = element_text(face="plain"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
pred_reptiles_plot
#####

# Get raster predictions #####  ####
# Fig 1D) Predicted Clustered phylo for Amphibians ####
prep= c("roughness","sand" ,"bio4","bio5","bio13","bio15")
all_anf_gdm_table = GDM_DF_function(turnover_matrix = All_an_phylo.beta.sim_Full,
                                    weights_mat1 = pd_sd_an,
                                    best_model = TRUE,
                                    var_to_final_model= prep)
gdm_all_anf <- gdm::gdm(all_anf_gdm_table, geo=TRUE)

preds=stack(elevation_c, slope_roughness_c,relief_roughness_c , precip_rast,temp_rast,soil_raster)
rastTrans_all_anf <- gdm::gdm.transform(gdm_all_anf, subset(preds,prep))
rastTrans_all_anf<- raster::mask(rastTrans_all_anf,cerrado)
plot(rastTrans_all_anf)

sampling_pts = rasterToPoints(rastTrans_all_anf[[1]])
sampling_pts=sampling_pts[,-3]
where_about = raster::extract(rastTrans_all_anf,sampling_pts)
sampling_pts_2 = cbind(sampling_pts,where_about)
sampling_pts_3 = na.exclude(sampling_pts_2)
attr(sampling_pts_3, "ATT") <- NULL
head(sampling_pts_3)

for_agnes = sampling_pts_3[,-c(1:2)]
tree_an = agnes(for_agnes, method = "average")
cluster_membership=cutree(tree_an, k = 12)
cluster_membership2=cutree(tree_an, k = 12)

pts_to_plot= as.data.frame(cbind(sampling_pts_3[,c(1:2)],cluster_membership))
pts_to_plot$cluster_membership2 = cluster_membership2
coordinates(pts_to_plot)=pts_to_plot[,c(1:2)]
pred_rast = raster::rasterize(pts_to_plot,rastTrans_end_rep,field=pts_to_plot@data$cluster_membership)
pred_pol = st_as_sf(rasterToPolygons(raster::rasterize(pts_to_plot,rastTrans_end_rep,field=pts_to_plot@data$cluster_membership2),dissolve = TRUE))
pts_to_plot_df = pts_to_plot@data
plot(pred_rast)
plot(pred_pol,fill=NA, add=TRUE)

## Getting a dendrogram ####
dend <- as.dendrogram(tree_an)

# Checking the order of the colnames ###
head(colnames(as.matrix(turnover_matrix)))
head(cluster_membership)
order.dendrogram(dend)

# Changing cluster membership according to the dend order
cluster_membership3=cluster_membership2[order(match(names(cluster_membership2), labels(dend)))]

dend2=branches_attr_by_clusters(dend,cluster_membership2)

# This is to visualize to which branches are the clusters assigned
dend4= dend2
labels(dend4) = cluster_membership2
dend4 %>% set("branches_lwd",2.5) %>% ladderize(TRUE)  %>%  plot(horiz  = FALSE,labels=FALSE,axes=FALSE)

# to get the right order from the phylogenetic turnover mapping
unique(cluster_membership2)
unique(get_leaves_branches_col(dend3))

# To map values in the raster (pts for ggplot)!
pts_to_plot_df$col= plyr::mapvalues(pts_to_plot_df$cluster_membership2, from=unique(cluster_membership3), to=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Color the dendrogram
dend3 <- color_branches(dend2,clusters = cluster_membership3,  col=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#5ab4ac','#abd9e9','#74add1','#4575b4','#313695',"darkblue"))

# Testing the order of the colors
#unique(pts_to_plot_df$cluster_membership)

# Produce a ggplot to be added on the above map!
par(bg="transparent")
dude_pred_an = dend3 %>% set("branches_lwd",.4) %>% set("labels", NA)  %>% ggplot(horiz  = TRUE,labels=FALSE,axes=FALSE)
#ggdend(dend3, horiz  = TRUE,labels=FALSE,axes=FALSE)

cols= unique(get_leaves_branches_col(dend3))
#cols=unique(pts_to_plot_df$col)


# Plotting itself #####
expression_title = paste0("(d) Amphibians- Predicted Phylogenetic Turnover")
pred_anf_plot= ggplot(map_data_sf)+
  geom_sf(data=map_data_sf,fill="grey60",
          col="grey30",
          size=0.1)+
  #geom_sf(data=cerrado_sf,colour="grey",fill="white", size=.25)+
  #geom_point(aes(x=x, y=y), colour = pts_df$col,data=pts_df,size=3.5)+
  geom_raster(data =pts_to_plot_df,
              aes(x=x,
                  y=y),
              fill= pts_to_plot_df$col)+
  geom_sf(data=map_data_sf,
          colour="grey30",
          fill=NA,size=0.1)+
  geom_sf(data=pred_pol,colour="black",fill=NA,size=0.2)+
  coord_sf(xlim=c(min_x,max_x),ylim = c(min_y,max_y))+
  scale_x_continuous(breaks = seq(-70, -40, by = 10))+
  scale_y_continuous(breaks = seq(-30, 10, by = 10))+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'white',size = .5),
        panel.background = element_rect(fill = 'grey95'))+
  theme(panel.border = element_rect(colour = NA, fill=NA, size=.5))+
  theme(axis.text=element_text(size=8))+
  theme(axis.line = element_line(colour = 'grey20', size = .1))+
  theme(legend.background=element_rect(fill = NA, colour = NA),
        legend.title = element_text(color="grey20",family="Times"),
        legend.text  = element_text(color="grey20",family="Times"),
        legend.key.height=unit(0.04,"npc"),
        legend.key.width=unit(0.01,"npc"),
        legend.position = c(1.0001,.75),
        legend.justification = 'left')+
  ggtitle(expression_title)+ 
  geom_rect(xmin = -Inf, xmax = -5050000,   ymin = 370000,    ymax = -1400000, fill = "grey60") +
  annotate(geom="text", x=-5800000, y= 295000, label="UPGMA clustering",
           color="white",size = 3.5,family="Times")+
  theme(plot.title = element_text(size=10,face="plain",hjust = 0,family="Times"))+
  theme(plot.title = element_text(face="plain"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
pred_anf_plot

## South America Insert
SA_MAP = ggplot(map_data_sf)+
  geom_sf(data=map_data_sf,fill="grey60",
          col="grey30",
          size=0.1)+
  geom_sf(data=cerrado_sf,colour="grey",fill="white", size=.25)+
  coord_sf(xlim=c(-8123756,max_x),ylim = c(-5952196,1552196))+
  theme_bw()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.background=element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        panel.background = element_rect(fill = 'white'))



#####
## Plotting hierarquical turnover ####
## Plotting side by side All ####
# Defining the positions and sizes of each map
vpa_ <- viewport(width = unit(8.3,"cm"), height = unit(9,"cm"), x = 0.25, y = 0.75)  
vpb_ <- viewport(width = unit(8.3,"cm"), height = unit(9,"cm"), x = 0.7, y = 0.75)  
vpc_ <- viewport(width = unit(8.3,"cm"), height = unit(9,"cm"), x = 0.25, y = 0.25)  
vpd_ <- viewport(width = unit(8.3,"cm"), height = unit(9,"cm"), x = 0.7, y = 0.25) 
vph_ <- viewport(width = unit(3.4,"cm"), height = unit(4.5,"cm"), x = 0.150, y = 0.83)
vph2_ <- viewport(width = unit(3.4,"cm"), height = unit(4.5,"cm"), x = 0.6, y = 0.83)
vph3_ <- viewport(width = unit(3.4,"cm"), height = unit(4.5,"cm"), x = 0.150, y = 0.33)
vph4_ <- viewport(width = unit(3.4,"cm"), height = unit(4.5,"cm"), x = 0.6, y = 0.33)
vpSA <- viewport(width = unit(2.8,"cm"), height = unit(3,"cm"), x = 0.865, y = 0.89)

# To save as a PDF file
#pdf(file=paste0("/Users/josue/Dropbox/1Doutorado/Chapter 1/cerrado_biogeo3/output_figures/",gsub(":", "-", Sys.time()),"Hclut_PhyBeta_imp.pdf"),onefile = TRUE)

# Open a new plotting device (deactivate this line if using the pdf function above)
quartz()
# FIG A
print(reptiles_plot ,vp = vpa_)
print(dude+
        theme(plot.background = element_rect(fill = "transparent", color = NA)), vp = vph_)

# FIG B
print(anf_plot, vp = vpb_)
print(dude_an+
        theme(plot.background = element_rect(fill = "transparent", color = NA)), vp = vph2_)

# FIG C
print(pred_reptiles_plot ,vp = vpc_)
print(dude_pred + theme(plot.background = element_rect(fill = "transparent", color = NA)) ,vp = vph3_)

# FIG D
print(pred_anf_plot ,vp = vpd_)
print(dude_pred_an + theme(plot.background = element_rect(fill = "transparent", color = NA)) ,vp = vph4_)

# South America insert
print(SA_MAP ,vp = vpSA)

dev.off()


#######
#End
