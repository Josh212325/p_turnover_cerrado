# Scripts to run GDM analyzes from Azevedo et al., 2021
# Here you you'll find how to produce Table 3 and Figure 1 and how to calculate total contributions of geography, environment and geo+env in GDM models 
# Author: Josué A. R. Azevedo
# Department of biological and environmental sciences - Uni of Gothenburg
# GGBC Sweden
# INPA AM Brazil
# email: josueanderson21@gmail.com

library(cluster)
library(ecodist)
library(robustbase)
library(plyr)
library(vegan)
library(betapart)
library(pbapply)
library(gdm)
library(parallel)
library(plyr)
library(rowr)
library(reshape)
library(ggplot2)
library(grid)

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

# Phylogenetic turnover was calculated with the function phylo.beta.pair from the R package betapart as follows:
# Phylo-beta-dissimilarity ###
# All reptiles all #
#hey <-function(sample_of_100_phylogenies){
#  b <-phylo.beta.pair(t(any_presence_abs_matrix),sample_of_100_phylogenies)
#  return(b$phylo.beta.sim)
#}
#hehey <-  pblapply(sample_of_100_phylogenies,hey,cl=7) 
#bla<- nrow(as.data.frame(as.matrix(hehey[[1]])))
#matricizando <- lapply(hehey,as.matrix)
#tmp <- plyr::ldply(matricizando) # convert to df
#tmp$counter <- 1:bla # 
#all.data1<- plyr::ddply(tmp, .(counter), function(x) colMedians(as.matrix(x[1:ncol(x)])))
#All_sq_phylo.beta.sim_FULLY <- as.dist(as.matrix(all.data1[,all.data1$counter]))

All_sq_phylo.beta.sim_FULLY = cerrado_layers[[2]]  # matrix of phylogenetic turnover for reptiles

Endemic_sq_phylo.beta.sim_FULLY = cerrado_layers[[3]] # matrix of phylogenetic turnover for endemic reptiles 

Endemic_an_phylo.beta.sim_Full = cerrado_layers[[4]]  # matrix of phylogenetic turnover for endemic amphibians

All_an_phylo.beta.sim_Full = cerrado_layers[[5]]  # matrix of phylogenetic turnover for amphibians

# Taxonomic turnover was calculated with the function beta.pair from the R package betapart:
#all_sq.betapair<-beta.pair(t(any_presence_abs_matrix), index.family="sor")
#all_sq.beta.sim <- as.dist(all_sq.betapair$beta.sim)

all_an.betapair = cerrado_layers[[6]]  # matrix of taxonomic turnover for amphibians

END_an.beta.sim = cerrado_layers[[7]]  # matrix of taxonomic turnover for endemic amphibians

all_sq.beta.sim = cerrado_layers[[8]] # matrix of taxonomic turnover for reptiles

END_sq.beta.sim = cerrado_layers[[9]]  # matrix of taxonomic turnover for endemic amphibians

# Standard deviation in phylogenetic diversity was calculated as follows
#hey <-function(sample_of_100_phylogenies){
#  pd_dist =picante::pd(t(any_presence_abs_matrix), sample_of_100_phylogenies, include.root=TRUE)
#  return(pd_dist$PD)
#}
#hehey <-  pblapply(sample_of_100_phylogenies,hey,cl=7) 
#bla<- nrow(as.data.frame(as.matrix(hehey[[1]])))
#matricizando <- lapply(hehey,as.matrix)
#tmp <- plyr::ldply(matricizando) # convert to df
#tmp$counter <- 1:bla # 
#all.data1<- plyr::ddply(tmp, .(counter), function(x) matrixStats::colSds(as.matrix(x[1:ncol(x)])))
#pd_sd_an = cbind(colnames(end_an_ten_or_more_an) ,as.data.frame(all.data1$V1),as.data.frame(pd_dist_test_an$SR))
#colnames(pd_sd_an)=c("site","pd_sd","SR")

pd_sd_sq = cerrado_layers[[10]] # standard deviation of PD per location across a sample of 100 phylogenies for reptiles

pd_sd_an = cerrado_layers[[11]] # standard deviation of PD per location across a sample of 100 phylogenies for amphibians
######

# Defining Behrmann equal area projection
behr <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

#####

# Running GDM #######
# Here I standardized the code by creating a function to run with any of the turnover matrices and environmental variables
#####

# Function parameters ####
turnover_matrix <- Endemic_sq_phylo.beta.sim_FULLY # turnover matrix
weights_mat1 <- pd_sd_sq # weighting matrix (ex., species richness, PD std)
best_model = FALSE # If true, it will select only the variables below
var_to_final_model= c("bio5","bio15") # Variables selected in the best model
#####

#### GDM_DF_function: a function to produce the required GDM tables.
#####

GDM_DF_function= function(turnover_matrix,weights_mat1,best_model,var_to_final_model){
  
  # Preparing the turnover matrix from a dist object 
  tax_tur <- as.data.frame(as.matrix(turnover_matrix))
  tax_tur=tax_tur[rownames(tax_tur) %in% grid_todos@data$layer,colnames(tax_tur) %in% grid_todos@data$layer]
  tax_tur<- cbind(row.names(tax_tur),tax_tur)
  colnames(tax_tur) <- c("site",  colnames(tax_tur[2:ncol(tax_tur)]))
  tax_tur[1:5, 1:5]
  
  # Getting the centroids of each grid cell
  grid_cells=grid_todos[grid_todos@data$layer %in% tax_tur$site,]
  coodinates_grid_cells <- rgeos::gCentroid(grid_cells, byid=TRUE,id =grid_cells@data$layer)
  coodinates_grid_cells <- as.data.frame(coodinates_grid_cells@coords)
  coodinates_grid_cells$pred <- rep(1 , nrow(coodinates_grid_cells))
  
  ## Preparing sites
  site <- tax_tur$site
  fake_preds <- cbind(row.names(tax_tur),as.data.frame(coodinates_grid_cells$pred))
  fake_preds$x <- coodinates_grid_cells$x
  fake_preds$y <- coodinates_grid_cells$y
  colnames(fake_preds) <- c("site", "fake_var","x","y")
  
  ## Adding weights (Standard deviation of PD)
  #weights_mat= weights_mat1[weights_mat1$site %in% tax_tur$site,]
  #weights_mat=weights_mat[,c("site","SR")]
  #colnames(weights_mat)=c("site","weights")
  
  # Preparing rasters of environmental variables
  preds=stack(elevation_c, slope_roughness_c,relief_roughness_c , precip_rast,temp_rast,soil_raster)
  preds= mask(preds,temp_rast[[1]])
  
  # Or if you wanted to use all CHELSA variables (pre-selected in Azevedo et al., 2021 using usdm::vifstep)
  #bio23=stack("/path/to/your/variables/CHELSA.tif")
  #names(bio23)=c("bio1","bio10", "bio11",'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9')
  #bio24= raster::aggregate(bio23,2,median) # Adjust to the resolution you want to use
  #bio25 <- crop(bio24,limits_models)
  #bio26 <- mask(bio25,limits_models)
  #bio27 <- resample(bio26,soil_raster[[1]])
  #preds = stack(elevation_c, slope_roughness_c,relief_roughness_c ,bio27,soil_raster)
  
  if (best_model==TRUE){
    preds = subset(preds, var_to_final_model)
  } else { print("full_model")}
  #plot(preds)
  names(preds)
  
  #preds=spatialEco::raster.transformation(preds, trans ="std") # If you wanted to standardize all variables
  preds_match=raster::extract(preds,grid_cells,fun="median")
  preds_match2= cbind(preds_match,fake_preds[,c(1,3:4)])
  #tail(preds_match2)
  #summary(preds_match2)
  
  # TABLE with weights (PD standard deviation)
  #gdmTab.dis <- formatsitepair(tax_tur, bioFormat=3, XColumn="x", YColumn="y",predData=preds_match2, siteColumn="site",weightType="custom",custWeights=weights_mat)
  
  # TABLE without weights
  gdmTab.dis <- formatsitepair(tax_tur, bioFormat=3, XColumn="x", YColumn="y",predData=preds_match2, siteColumn="site")
  return(gdmTab.dis)
}
#####

# Listing all phylogentic turnover tables in a single object (for running in a loop)
all_diss_talbes= list(All_sq_phylo.beta.sim_FULLY,
                      Endemic_sq_phylo.beta.sim_FULLY,
                      Endemic_an_phylo.beta.sim_Full,
                      All_an_phylo.beta.sim_Full)

length(all_diss_talbes)

# And here I will run the variable importance function of the gdm R package (gdm.varImp) for all phylogenetic turnover matrices. This may take one entire day (or more) to run with this dataset. If only testing this script, I would recommend to run only 10 permutations if you do not have a powerful computer  ####
# If using weights_mat1, create a list of of weighting matrices - list(pd_sd_sq,pd_sd_sq,pd_sd_an,pd_sd_an) - in the same other of your phylogenetic turnover list (all_diss_talbes), activate the option for weights in the GDM_DF_function above, add this list to the function arguments below, and use mapply instead of lapply ####

best_model=FALSE # For running all possible combination of variables
# Function to run variable selection for all P_turnover tables
getting_all_diss_models=function(all_diss_talbes){
  gdmTab.dis=GDM_DF_function(turnover_matrix=all_diss_talbes,best_model=FALSE)
  var_importance= gdm.varImp(gdmTab.dis, geo=TRUE, nPerm=10, parallel=TRUE, cores=5,sampleSitePairs=0.8)
  return(var_importance)
}
model_selection_all1000 = lapply(all_diss_talbes,getting_all_diss_models)
#####

#
# If you do not run the function gdm.varImp, go directly to the section: "Here is after model selection" ####

#####

# Selecting the model with the highest value of explained deviance which retained only variables that were important in more than 160 rounds of permutation ####

i=4
df_model_sel_all= c(NULL)
for (i in c(1:4)){
  x = as.data.frame(model_selection_all1000[[i]][[3]] < 0.160)
  
  df_model_sel = head(colnames(x[sapply(x, function(x) !any(which(!x)))]),1)
  
  df_model_sel_all = c(df_model_sel_all,df_model_sel)
}

df_model_sel_all= cbind(df_model_sel_all,c("All_sq_phylo.beta.sim_FULLY",
                                           "Endemic_sq_phylo.beta.sim_FULLY",
                                           "Endemic_an_phylo.beta.sim_Full",
                                           "All_an_phylo.beta.sim_Full"))

colnames(df_model_sel_all)= c("model","dataset")
df_model_sel_all=as.data.frame(df_model_sel_all)


# Table 3 (Azevedo et al., 2021), first part (Best Models = Model deviation, Percentage of deviation explained) ####
i=1
which_group = c("Reptiles", "Endemic Reptiles","Endemic Amphibians","Amphibians")
part_a= NA # 
for (i in c(1:4)){
  bm= df_model_sel_all[i,1]
  model_x = round(as.data.frame(model_selection_all1000[[i]][[1]][,bm]),1)
  colnames(model_x) = which_group[i]
  part_a= cbind(part_a,model_x)
}
part_a=part_a[c(1,2,4),c(2,3,5,4)]
write.csv(part_a) #You may want to copy and past the content of part_a displayed in the console in the Microsoft Word and directly transform the results on a table

######

# For the second part of Table 3 (Azevedo et al 2021). Getting Variable importance ####

model_y_rows = c(NULL) # run this first one time
for (i in c(1:4)){ # 
  bm= df_model_sel_all[i,1]
  model_y = round(as.data.frame(model_selection_all1000[[i]][[2]][,paste0(bm)]),1)
  model_y$sig= as.data.frame(model_selection_all1000[[i]][[3]][,paste0(bm)])
  #model_y$sig= ifelse(model_y$sig < 0.05,"",NA )
  model_y = na.exclude(model_y)
  model_y=rownames(model_y)
  model_y_rows = c(model_y_rows,model_y)
}
model_y_rows = unique(model_y_rows)
initial_df = as.data.frame(cbind(model_y_rows,rep(1,length(model_y_rows))))
colnames(initial_df)= c("variable","dummy_var")

part_b= initial_df # run this first one time
for (i in c(1:4)){ # 
  bm= df_model_sel_all[i,1]
  model_y = round(as.data.frame(model_selection_all1000[[i]][[2]][,paste0(bm)]),1)
  model_y$sig= as.data.frame(model_selection_all1000[[i]][[3]][,paste0(bm)])
  #model_y$sig= ifelse(model_y$sig < 0.05,"",NA )
  model_y = na.exclude(model_y)
  model_y$var =rep(rownames(model_y))
  colnames(model_y) = c("Variable importance","sigui","variable")
  attr(model_y[,2], "ATT") <- NULL
  p_value_sar_res= ifelse(model_y$sigui >= 0.05, p_value_sar_res <- paste0(""), ifelse(model_y$sigui <= 0.05 & model_y$sigui >= 0.01, p_value_sar_res <- paste0("*"), ifelse(model_y$sigui < 0.009 & model_y$sigui >= 0.001, p_value_sar_res <- paste0("**"),  p_value_sar_res <- paste0("***")))) ## p_value_sar refers to the asterisks 
  model_y$`Variable importance`= paste0(model_y$`Variable importance`,p_value_sar_res)
  model_y=model_y[,c("variable","Variable importance")]
  colnames(model_y) = c("variable",paste0(which_group[i]))
  part_b = merge(part_b,model_y, by="variable",all.x=TRUE)
}
part_b= part_b[c(6,7,8,4,5,1,2,3),c(1,3,4,6,5)] # reordering rows
# reordering DF
write.csv(part_b) #You may want to copy and past the content of part_b displayed in the console in the Microsoft Word and directly transform the results on a table


#############################################
##### Here is after model selection #########
#############################################

# Assigning the variables selected in the best model for each group the the next steps
var_to_final_model_rep = c("bio5","bio15")
modTest_anf = c("roughness","sand" ,"bio4","bio5","bio13","bio15")
var_to_final_model_end_rep = c("bio13","bio15","bio18")
modTest_en_anf = c("roughness","sand" ,"bio5","bio13","bio15")


#####

###################################################
## Figure 1 #######################################
###################################################

# Here I created a function (extract_splines) for getting splines for each variable. Note that lines in all_DF_x are particular to my dataset  #####
extract_splines= function(gdm_table,group_name){
  model= gdm::gdm(gdm_table, geo = TRUE)
  gdm_rep_splineDat= gdm::isplineExtract(model)
  all_DF_y = cbind.data.frame(rownames(as.data.frame(gdm_rep_splineDat$y)) ,gdm_rep_splineDat$y)
  colnames(all_DF_y) = c("axis",colnames(all_DF_y)[-1])
  all_DF_y =melt(all_DF_y, variable.name= "axis" )
  colnames(all_DF_y) = c("axis", "varible","y_axis")
  
  all_DF_x = cbind.data.frame(rownames(as.data.frame(gdm_rep_splineDat$x)) ,gdm_rep_splineDat$x)
  colnames(all_DF_x) = c("axis",colnames(all_DF_x)[-1])
  all_DF_x$Geographic=all_DF_x$Geographic/1000
  all_DF_x$bio4 = all_DF_x$bio4/1000
  all_DF_x$bio5 = all_DF_x$bio5/10
  all_DF_x$bio15 = all_DF_x$bio15/10
  #all_DF_x$bio15 = all_DF_x$bio11/10
  all_DF_x =melt(all_DF_x, variable.name= "axis" )
  colnames(all_DF_x) = c("axis", "varible","x_axis")
  
  all_DF_2 = cbind.data.frame(all_DF_y[,-1],all_DF_x[,-c(1:2)])
  colnames(all_DF_2) = c("varible", "y_axis","x_axis")
  
  test=all_DF_2[all_DF_2$varible=="sand",]
  
  tail(all_DF_2)
  all_DF_3=cbind(all_DF_2, rep(paste0(group_name),nrow(all_DF_2)))
  head(all_DF_3,20)
  colnames(all_DF_3) = c("variable","y_axis","x_axis", "order")
  tail(all_DF_3)
  return(all_DF_3)
}

# Most important variables selected in all four best models.
var_to_final_model = c("roughness","sand","bio4","bio5","bio13","bio15","bio18")

#"bio4"=Temperature Seasonality
#"bio5"=Max Temperature of Warmest Month
#"bio13"=Precipitation of Wettest Month
#"bio15"=Precipitation Seasonality
#"bio18"=Precipitation of Warmest Quarter

# Run GDM_DF_function to include the above variables for each Order and extract the splines

# Reptiles
rep_best_gdm_table = GDM_DF_function(turnover_matrix = All_sq_phylo.beta.sim_FULLY,
                                     weights_mat1 = pd_sd_sq,
                                     best_model = TRUE,
                                     var_to_final_model= var_to_final_model)

rep_splines= extract_splines(gdm_table=rep_best_gdm_table,group_name="Reptiles")
tail(rep_splines)


# Endemic Reptiles
End_rep_best_gdm_table = GDM_DF_function(turnover_matrix = Endemic_sq_phylo.beta.sim_FULLY,
                                         weights_mat1 = pd_sd_sq,
                                         best_model = TRUE,
                                         var_to_final_model= var_to_final_model)

end_rep_splines= extract_splines(gdm_table=End_rep_best_gdm_table,group_name="Endemic_Reptiles")
tail(end_rep_splines)

# Amphibians
anf_best_gdm_table = GDM_DF_function(turnover_matrix = All_an_phylo.beta.sim_Full,
                                     weights_mat1 = pd_sd_an,
                                     best_model = TRUE,
                                     var_to_final_model= var_to_final_model)
Amph_splines= extract_splines(gdm_table=anf_best_gdm_table,group_name="Amphibians")
tail(Amph_splines)

# Endemic Amphibians
End_anf_best_gdm_table = GDM_DF_function(turnover_matrix =Endemic_an_phylo.beta.sim_Full,
                                         weights_mat1 = pd_sd_an,
                                         best_model = TRUE,
                                         var_to_final_model= var_to_final_model)
End_Amph_splines= extract_splines(gdm_table=End_anf_best_gdm_table,group_name="Endemic Amphibians")
tail(End_Amph_splines)
#####

# Merging all dataframes of xy spline values    ####
all_splines = rbind(rep_splines,end_rep_splines,Amph_splines,End_Amph_splines)

head(all_splines)
tail(all_splines)
str(all_splines)

unique(all_splines$variable) # Always check the order of your variables here to set the name of the variables in the next step

#"bio4"=Temperature Seasonality
#"bio5"=Max Temperature of Warmest Month
#"bio13"=Precipitation of Wettest Month
#"bio15"=Precipitation Seasonality
#"bio18"=Precipitation of Warmest Quarter

all_splines$variable = plyr::mapvalues(all_splines$variable, from=unique(all_splines$variable), to=c("Geographic (km)","Relief Roughness (m)","Sand (%)","Temperature \n Seasonality (°C)","*Max Temperature \n Warmest Month (°C)", "Precipitation \n Wettest Month (mm)", "Precipitation \n Seasonality (mm)","**Precipitation \n Warmest Quarter (mm)"))

unique(all_splines$variable) # Always check this order

#all_splines = all_splines[order(factor(all_splines$variable,levels=c(c("Geographic","Elevation","Relief Roughness","Sand","Temperature \n Seasonality","Max. Temp. \n Warmest Month","Precipitation of \n Wettest Month","**Mean Temp. \n Coldest Quarter","Precipitation \n Seasonality")))),]

#####

# Plotting All splines, Fig 1 Azevedo et al., 2021 ####
quartz()

# For all spp
bp <- ggplot(all_splines, aes(x=x_axis, y=y_axis, group=order))+
  geom_line(aes(colour=order),size=1.25)+
  scale_color_manual(values = c('#d7191c','#fdae61','#2c7bb6','#abd9e9'),labels=c("Reptiles","Endemic Rept.","Amphibians","Endemic Amph."),name="Orders")+
  #scale_color_manual(values = c('#d7191c','#2c7bb6'),labels=c("Reptiles","Amphibians"),name="Orders")+
  #scale_color_manual(values = c('#d7191c','#fdae61','#2c7bb6','#abd9e9'),labels=c("Endemic Rept.","Endemic Amph."),name="Orders")+
  ylab("Phylogenetic Turnover (Partial ecological distance)")+
  xlab(NULL)+
  theme_light()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,hjust = 1),
        text = element_text(family="Times",size = 13))+
  facet_wrap(variable ~ .,scales="free_x",ncol = 4,strip.position="bottom")+
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=10.5,colour = "black",vjust=1.2),strip.placement = "outside")
bp

# To save as a PFD file
pdf(file=paste0("/Users/josue/Dropbox/1Doutorado/Chapter 1/cerrado_biogeo3/output_figures/",gsub(":", "-", Sys.time()),"All_splines_2.pdf"),onefile = TRUE,width=6.61417,height=9.0551181) ##
#168mm width per 230mm of height
vpa_ <- viewport(width = unit(16.8,"cm"), height = unit(18,"cm"), x = 0.5, y = 0.5)
print(bp,vp=vpa_)
dev.off()


####################################
# Calculating model contribution (only geography, only enviroment, geo+env) ##
####################################

# Function to extract deviance explained & calculate unique & shared contributions
# Look that this step is particular to the way I named my variables in -grep
extract_importance = function(gdm_table,group_name){
  #Movel all
  model_all= gdm::gdm(gdm_table, geo = TRUE)
  #Model env
  model_Env= gdm::gdm(gdm_table, geo = FALSE)
  #Movel space only
  gdm_table_XY = gdm_table[,-grep(pattern ="sand|roughness|bio",colnames(gdm_table),value=FALSE)]
  colnames(gdm_table_XY)
  model_XY= gdm::gdm(gdm_table_XY, geo = TRUE)
  
  model_all$explained
  model_Env$explained
  model_XY$explained
  
  gdm..devAll <- model_all$explained
  gdm..devEnv <- model_Env$explained
  gdm..devXY <- model_XY$explained
  
  # Calculations as in Fitzpatrick et al., 2013  
  gdm..pDevEnv <- (gdm..devAll-gdm..devXY)/gdm..devAll
  gdm..pDevXY <- (gdm..devAll-gdm..devEnv)/gdm..devAll
  gdm..pDevAll <- 1 - gdm..pDevEnv - gdm..pDevXY
  
  gdm..pDevAll+gdm..pDevXY+gdm..pDevEnv
  
  explayned_by = c("purely environment","purely space", "shared env & space") 
  
  values_prop = c(gdm..pDevEnv,gdm..pDevXY,gdm..pDevAll)
  
  proportions =  cbind.data.frame(explayned_by,round(values_prop,2),rep(group_name,3))
  
  return(proportions)
}

var_to_final_model_rep = c("bio5","bio15")
modTest_anf = c("roughness","sand" ,"bio4","bio5","bio13","bio15")
var_to_final_model_end_rep = c("bio13","bio15","bio18")
modTest_en_anf = c("roughness","sand" ,"bio5","bio13","bio15")

# Reptiles
rep_best_gdm_table = GDM_DF_function(turnover_matrix = All_sq_phylo.beta.sim_FULLY,
                                     weights_mat1 = pd_sd_sq,
                                     best_model = TRUE,
                                     var_to_final_model= var_to_final_model_rep)

var_exp_rep = extract_importance(rep_best_gdm_table,"Reptiles")
var_exp_rep

# Endemic Reptiles
End_rep_best_gdm_table = GDM_DF_function(turnover_matrix = Endemic_sq_phylo.beta.sim_FULLY,
                                         weights_mat1 = pd_sd_sq,
                                         best_model = TRUE,
                                         var_to_final_model= var_to_final_model_end_rep)
var_exp_end_rep = extract_importance(End_rep_best_gdm_table,"Endemic Rep.")

# Amphibians
anf_best_gdm_table = GDM_DF_function(turnover_matrix = All_an_phylo.beta.sim_Full,
                                     weights_mat1 = pd_sd_an,
                                     best_model = TRUE,
                                     var_to_final_model= modTest_anf)

var_exp_anf = extract_importance(anf_best_gdm_table,"Amphibians")


# Endemic Amphibians
End_anf_best_gdm_table = GDM_DF_function(turnover_matrix =Endemic_an_phylo.beta.sim_Full,
                                         weights_mat1 = pd_sd_an,
                                         best_model = TRUE,
                                         var_to_final_model= modTest_en_anf)

var_exp_end_anf = extract_importance(End_anf_best_gdm_table,"Endemic amphib.")

var_exp_all = rbind.data.frame(var_exp_rep,var_exp_end_rep,var_exp_anf,var_exp_end_anf)
colnames(var_exp_all)=c("explayned_by","Perc","group")


var_exp_all <- ddply(var_exp_all, .(group,Perc), 
                     transform, pos = cumsum(Perc) - (0.5 * Perc))

# Create the barplot (Figure in supplementary material -Azevedo et al., 2021)
library(ggthemes)
ggplot(data=var_exp_all, aes(x=group, y=Perc,label = Perc, fill=explayned_by)) +
  geom_bar(stat="identity")+
  geom_text(size = 3,col="white", position = position_stack(vjust = 0.5))+
  #scale_fill_brewer(palette="Accent")+
  scale_fill_manual(values=c('#d6604d','#4393c3','#b2abd2'),name="Explained by")+
  ylab("Proportion of total exaplained deviance")+
  xlab(NULL)+
  theme_tufte()

# End
