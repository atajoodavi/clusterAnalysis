# Title: cluster.R
# Author: Nathan Johnson
# Date: 2014.06.20

###  REQUIRED PACKAGES  ###
require(reshape2)
require(gdata) # nobs() function
require(cluster)
require(clusterSim)
require(zoo)
require(maptools)
require(maps)
require(mapdata)
require(mapproj)
require(lubridate) # year() function
library(pvclust)
# detach(package:mclust)

rm(gap) # no plot of gaps from old data

mainDir = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/"
new.data = "Y" # "Y", "N"
filled.Data = "Filled" ## "Filled", "Original"

if(new.data == "Y"){
  if(filled.Data =="Original"){
    model.domain.data1 = read.csv(paste("C:/all.agency.mdata.csv", sep = ""), header = TRUE) ### all agency mdata
  } else {model.domain.data1 = read.csv(paste(mainDir, "data/Filled_stations_stack_second_All_agency_MaxCorr_0.9_MatchPair_10_MinRegPeriod3.csv", sep = ""), header = TRUE)} ### all agency mdata
}
model.domain.data = model.domain.data1

## USER INPUT ##
model.area = "GA" # which region for cluster analysis ("CSM", "NFSEG", "All_agency", "GA")
aquifer.type = c("UFA") # "UFA", "MultiAquifer", "SAS", "ICU", "UZLFA", "check", "LFA", "ULFA", "FPZ", "MCU", "Crystalline Ride Aquifer", "Brunswick Aquifer System", "LSCU", "SECPA", "Bottom Aquf", "Crystalline Rock Aquifer", "Valley and Ridge Aquifer"
GIS.inventory = read.csv(paste(mainDir, "data/inventory/All_agency_well_inventory_strat.csv", sep =""))

if(model.area != "All_agency"){
  GIS.inventory = subset(GIS.inventory, !is.na(GIS.inventory[model.area]))}

model.domain.data = subset(
  model.domain.data, subset = model.domain.data$univName %in% GIS.inventory$univName)

if(aquifer.type != ""){ # create column with the aquifer type
  GIS.inventory.aquifer = subset(
    GIS.inventory, GIS.inventory$aquiferFinal %in% aquifer.type) # subset stations based on aquifer.type
  model.domain.data = subset(model.domain.data, model.domain.data$univName %in% GIS.inventory.aquifer$univName)
}

#### DEFINED PARAMETERS ####
# station.list = c("SJRWMD15282844","SJRWMD19324384") # Well selected for greatest correlation with another
station.var = "univName"
level.var = "levelNAVD88"
date.var = "monthYear"
start.date <- as.Date("01/01/1982",format = "%m/%d/%Y")  # Begin date for the analysis
end.date <- as.Date("12/31/2012",format = "%m/%d/%Y")  # End data for the analysis
percent.threshold = 0.9 # percent of data available over the period selected
cluster.number = 21 #number of clusters to design

if(model.area == "CSM"){
  map.xlim = c(-84,-80.5) # CSM longitute min and max
  map.ylim = c(27.5,30) # CSM latitude min and max
  county.lines = "Y" # display county lines on map ("Y" or "N")
  png.width = 650
  png.height = 550
} else if(model.area == "NFSEG"){
  map.xlim = c(-85,-80) # NFSEG longitute min and max
  map.ylim = c(28.5,33.5) # NFSEG latitude min and max
  county.lines = "N" # display county lines on map ("Y" or "N")
  png.width = 650
  png.height = 700
} else if(model.area == "All_agency"){
  map.xlim = c(-86,-80) # All agency longitute min and max
  map.ylim = c(25,34) # All agency latitude min and max
  county.lines = "N" # display county lines on map ("Y" or "N")
  png.width = 650
  png.height = 1000
} else if(model.area == "GA"){
  map.xlim = c(-85.5,-80.5) # NFSEG longitute min and max
  map.ylim = c(30.3,33.5) # NFSEG latitude min and max
  county.lines = "N" # display county lines on map ("Y" or "N")
  png.width = 650
  png.height = 550
} else stop("Select appropriate model.area that has and inventory .csv file 
            with appropriate latDD and lonDD columns or create a 
            model domain in the previous if else statements")

####   DATA FORMATTING AND PERIOD SELECTION   ####
wide.data = dcast(model.domain.data, monthYear~univName, value.var = "levelNAVD88")
wide.data$monthYear = as.Date(wide.data$monthYear, format = "%Y-%m-%d") # create data for monthYear
# wide.data$monthYear = as.Date(wide.data$monthYear, format = "%m/%d/%Y") # create data for monthYear
wide.data = wide.data[order(wide.data$monthYear),] # order the dataframe
wide.data = wide.data[which((wide.data$monthYear >= start.date) & 
                              (wide.data$monthYear <= end.date)),] # period of record 
# plot(wide.data$monthYear, wide.data[,8], main = "NWFWMD10121", ylab = "Level NAVD88", xlab = "Date")

####  SELECTION OF WELLS THAT HAVE DATA THAT MEET THRESHOLDS OF % COMPLETE AND PERIOD OF RECORD ####
col.good = numeric(0) # placeholder for numeric vector to record the stations that have adequate data

for (w in 2:length(wide.data[,2:length(wide.data[1,])])){
  num.row = length(wide.data[,w]) # total number of rows
  num.row.obs = nobs(wide.data[,w]) # counts the number of !NA values
  ratio.obs = (num.row.obs/num.row)
  if(ratio.obs > percent.threshold){
    col.good = c(col.good, w) # columns that meet the threshold
  }
}

if(length(col.good) == 1)(
  stop("There are no variables that meet this threshold"))
good.stations = data.frame(wide.data$monthYear,wide.data[,col.good]) # dataframe of stations that meet the threshold of % missing
colnames(good.stations)[1] =paste("date")
row.names(good.stations) = wide.data$monthYear

#### CREATE ZOO OBJECT TO FILL IN TIME SERIES GAPS ####
zoo.s = zoo(good.stations[2:length(good.stations[1,])],good.stations$date)
filled.monthly = na.approx(zoo.s, maxgap = 7, na.rm = FALSE) # fill in data with linear interpolation
filled.monthly = na.locf(filled.monthly, maxgap = 7, fromLast = TRUE) # fill in data with nearest real value from back to forward to fill values at beggining
filled.monthly = na.locf(filled.monthly, maxgap = 7, fromLast = FALSE) # fill in data with nearest real value forward to fill in end values
df.filled = data.frame(filled.monthly) # convert back to dataframe from zoo
df.filled = data.frame(data.frame(row.names(df.filled)), df.filled) # add date column to df.filled
colnames(df.filled)[1] = c("date")
df.filled = df.filled[order(df.filled$date),] # Sort df.filled by date

################### CLUSTER ANALYSIS ##################################
norm.df.filled = scale(df.filled[,2:length(df.filled[1,])]) #normalize the data by column
norm.df.filled = norm.df.filled[,colSums(is.na(norm.df.filled)) < nrow(norm.df.filled)] # remove columns where all values are NA
horz.norm.df.filled = as.data.frame(t(norm.df.filled)) # transpose to fit into function (agnes)
cluster = agnes(horz.norm.df.filled, diss = FALSE , stand = FALSE) #cluster analysis with good.stations

## CLUSTER LABELS FOR SPECIFIC NUMBER OF CLUSTERS ##
class = cutree(cluster, k = cluster.number)
df.class = data.frame(class,colnames(norm.df.filled)[1:length(colnames(norm.df.filled))])
colnames(df.class) = c("cluster_num", "univName")

## CREATE .CSV WITH CLUSTER NUMBERS ##
GIS.inventory.cluster = merge(x = GIS.inventory,y = df.class, by.x = "univName", by.y = "univName")
write.csv(GIS.inventory.cluster, 
          paste(mainDir, "data/clusters/",model.area,
                "_wells_clusters_",start.date,"-",end.date,"_",
                as.Date(Sys.time()),".csv", sep =""))

###### CREATE FIGURES ##########
subDir = paste("figures/cluster/", model.area, "/", sep = "")
png(paste(file = mainDir,subDir,"Cluster_plots_", filled.Data,"_",
          percent.threshold,"_", model.area,"_", aquifer.type, "_", 
          year(start.date), "-", year(end.date),"_Clus_num_", cluster.number,"_Run_", 
          as.Date(Sys.time()),".png", sep =""),
    width = 1000, height = 1000, units= "px", res=100)
par(mfrow = c(ceiling(sqrt(cluster.number+1)),ceiling(cluster.number/ceiling(sqrt(cluster.number+1)))))
prin.comp.all = data.frame(date = wide.data$monthYear)
prop.var = numeric()

## Cluster optimization gap figure ##
na.omit.df.filled = df.filled[,complete.cases(t(df.filled))] # remove stations with missing values
gap = clusGap(na.omit.df.filled[2:ncol(na.omit.df.filled)], FUNcluster = pam, K.max = 20, B= 3)
plot(gap, xlab = "Total Cluster Number", main = "Cluster Gap")

## Cluster figures ##
for(clus.number in 1:cluster.number){
  cluster.subset = subset(GIS.inventory.cluster, 
                          (GIS.inventory.cluster$cluster_num == clus.number))
  cluster.wells = data.frame(date = as.Date(df.filled[,1]), 
                             scale(df.filled[as.character(cluster.subset$univName)]),
                             stringsAsFactors= FALSE)

  plot(cluster.wells$date, cluster.wells[,2], type = "l", 
       xlab = "Date", ylab = "Stnd_Level", 
       main = (if(nrow(cluster.subset) == 1){paste("Cluster:",clus.number,aquifer.type, "\n", as.character(cluster.subset$univName), sep = " ")} else {paste("Cluster:", clus.number, aquifer.type)}), 
       col = 1, ylim = c(min(cluster.wells[,seq(2,(ncol(cluster.wells)))], na.rm = TRUE),
                         max(cluster.wells[,seq(2,(ncol(cluster.wells)))], na.rm = TRUE)))
  
  if(nrow(cluster.subset) > 1){
    for(p in 2:(nrow(cluster.subset))){
      col.name = as.character(cluster.subset[p,2])
      points(cluster.wells$date, cluster.wells[,(p+1)], 
             pch = p, col = p, type = "l")
    }}
  
  ### PRINCIPAL COMPONENT ANALYSIS ###
  na.omit.cluster.wells = cluster.wells[,complete.cases(t(cluster.wells))] # remove stations with missing values
  if(is.Date(na.omit.cluster.wells) == TRUE) {next} # when all stations contain Null values skip PCA
  scaled.wells = data.frame(scale(na.omit.cluster.wells[,2:ncol(na.omit.cluster.wells)], center = TRUE))
  
  if(ncol(scaled.wells) == 1){
    prin.comp.all = merge(prin.comp.all, na.omit.cluster.wells, by = "date", all = TRUE)
    prop.var = rbind(prop.var, 1)} else {
      prin.comp = prcomp(scaled.wells[1:ncol(scaled.wells)])
      if(prin.comp$rotation[1,1] < 0)(rotation = -1) else(rotation = 1) # correction for timeseries multiplied by sign of the rotation
      prin.comp1 = data.frame(date = as.Date(row.names(prin.comp$x)), 
                              prin.comp.ts = (rotation * scale(prin.comp$x[,1])))
      names(prin.comp1)[2] = paste("PC_", clus.number, sep = "")
      prin.comp.all = merge(prin.comp.all, prin.comp1, by = "date", all = TRUE)
      prop.var = rbind(prop.var, summary(prin.comp)$importance[2,1])
    }
  row.names(prin.comp.all) = wide.data$monthYear
}
dev.off()
if(prin.comp$rotation[1,1] < 0)(rotation = -1) else(rotation = 1)

### Cluster Analysis Map Figures ###
png(paste(file = mainDir,subDir,"Cluster_map_",filled.Data,"_",
          percent.threshold,"_",model.area,"_",aquifer.type, "_",  
          year(start.date), "-", year(end.date),"_Clus_num_", 
          cluster.number,"_Run_",as.Date(Sys.time()),".png", sep =""),
    width = png.width, height = png.height, units= "px", res=100)        
par(mfrow = c(1,1), mar = c(1,0,0,0))
if (county.lines == "Y"){map("county", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)
}else{map("state", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)}
map.axes()

leg.names.list = character()
for(clus.num in 1:cluster.number){
  cluster.subset = subset(GIS.inventory.cluster, GIS.inventory.cluster$cluster_num  == clus.num)
  points(x = cluster.subset$lonDD, y = cluster.subset$latDD, col = clus.num, pch = clus.num, cex = 0.75)
  leg.names = c(paste("Cluster:", clus.num))
  leg.names.list = c(leg.names.list, leg.names)}
legend("topleft", title = paste(aquifer.type," ",year(start.date), "-", year(end.date), sep = ""), 
       legend = leg.names.list, col = seq(1:length(leg.names.list)), 
       pch = seq(1:length(leg.names.list)), merge = FALSE, border=F, cex = 0.75)
dev.off()

### Principal Component Figures ###
subDir = paste("figures/principalComponent/", model.area, "/", sep = "")
png(paste(file = mainDir,subDir,"prinComp_plots_", filled.Data,"_",
          percent.threshold,"_", model.area,"_", aquifer.type, "_", 
          year(start.date), "-", year(end.date),"_Clus_num_",
          cluster.number,"_Run_", as.Date(Sys.time()),".png", sep =""),
    width = 1000, height = 1000, units= "px", res=100)
par(mfrow = c(ceiling(sqrt(cluster.number+1)),ceiling(cluster.number/ceiling(sqrt(cluster.number+1)))))
plot(gap, xlab = "Total Cluster Number", main = "Cluster Gap")
for(ploty in 2:ncol(prin.comp.all)){
  plot(prin.comp.all$date, prin.comp.all[,ploty], 
       main = paste("PC1 Cluster", (ploty-1), aquifer.type, sep = " "), type = "l",
       ylab = "Stnd_Comp",xlab = "Date") 
  legend ("topleft", legend = paste("Prop Var=", prop.var[(ploty-1),]), bty = "n")}
dev.off()

### PLOT WELLS IN INDIVIDUAL AQUIFER BY YEAR###
# subDir = paste("figures/aquiferStationMap/", model.area, "/",sep = "")
# aquifer.type = "UFA"
# png(paste(file = mainDir,subDir,model.area,"_",aquifer.type,"_byYear_Run_",as.Date(Sys.time()),".png", sep =""),
#     width = 650, height = 1000, units= "px", res=90,  antialias = "cleartype")
# par(mfrow = c(4,3), mar = c(1,0,0,0))
# for (year.type in 2000:2010){
# #     year.type = "2009"
#   data.year = paste("wl", year.type, sep = "")
#   if (county.lines == "Y"){map("county", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)
#   } else {map("state", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)}
#   map.axes()
#   aquifer.subset = subset(GIS.inventory, GIS.inventory$aquifer == aquifer.type)
#   aquifer.subset = subset(aquifer.subset, aquifer.subset[data.year] > 0)
#   points(x = aquifer.subset$lonDD, y = aquifer.subset$latDD, col = 1, pch = 3)
#   legend("bottomleft", legend = c(aquifer.type, year.type), merge = FALSE, border=F, cex = 1.2)
# }
# dev.off()

### PLOT WELLS IN INDIVIDUAL AQUIFER ###
# subDir = paste("figures/aquiferStationMap/", model.area, "/",sep = "")
# aquifer.type = c("SAS","UFA","ICU")
# all.data = "N" # if desired well coverage for whole period of record "Y" and "N" for wells selected that meet thresholds
# if(all.data == "N")(period = paste(" ",year(start.date), "-", year(end.date))) else (period = "POR")
# for(aq in 1:length(aquifer.type)){
#   png(paste(file = mainDir,subDir,model.area,"_",aquifer.type[aq],"_",filled.Data,"_",period,"_Run_",as.Date(Sys.time()),".png", sep =""),
#       width = 650, height = 1000, units= "px", res=90,  antialias = "cleartype")
#   if (county.lines == "Y"){map("county", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)
#   } else {map("state", xlim = map.xlim, ylim = map.ylim , col ="gray90", fill = TRUE)}
#   map.axes()
#   GIS.good.inventory = GIS.inventory[GIS.inventory$univName %in% colnames(good.stations), ]
#   if(all.data == "N")(
#     aquifer.subset = subset(GIS.good.inventory, GIS.good.inventory$aquiferFinal == aquifer.type[aq])
#   )else(aquifer.subset = subset(GIS.inventory, GIS.inventory$aquiferFinal == aquifer.type[aq]))
#   points(x = aquifer.subset$lonDD, y = aquifer.subset$latDD, col = 1, pch = 1)
#   legend("bottomleft", legend = paste(aquifer.type[aq],period, sep = " "), merge = FALSE, border=F, cex = 1.5)
#   dev.off()
# }

## SINGLE WELL MAP ##
# single.well = "SJRWMD11211584"
# map("state", xlim = map.xlim, ylim = map.ylim, col ="gray90", 
#     fill = TRUE, xlab = "Longitude (NAD83)", ylab = "Lattitude (NAD83)")
# map.axes()
# leg.names.list = character()
# cluster.subset = subset(GIS.inventory, GIS.inventory$univName  == single.well)
#   points(x = cluster.subset$lonDD, y = cluster.subset$latDD, col = "red", pch = 1)
#   leg.names = c(paste("Station:", single.well, sep = ""))
# legend("bottomleft", legend = leg.names, col = "red", 
#        pch = seq(1:length(leg.names.list)), merge = FALSE, border=F, cex = 0.75)

## Compare plots of orginal with filled ##
# station.select = 70
# plot.zoo(zoo.s[,station.select],
#       main = paste("Original ",colnames(filled.monthly)[station.select], sep = ""), 
#       ylab = "Level NAVD88 ft", xlab = "Date") 
# plot.zoo(filled.monthly[,station.select],
#           main = paste("Filled ",colnames(filled.monthly)[station.select],sep = ""), 
#           ylab = "Level NAVD88 ft", xlab = "Date")


## Cluster analysis Dendrogram
# library(pvclust)
# fit <- pvclust(df.filled[2:ncol(df.filled)], method.hclust="complete",
#    method.dist="correlation", use.cor = "complete.obs", nboot = 1000, r = c(0.5, 1, 1.5))
# plot(fit)
# pvrect(fit, alpha=.95)
# 
# fit1 = pvclust(na.omit.cluster.wells[2:ncol(na.omit.cluster.wells)],method.hclust="complete",
#    method.dist="correlation", use.cor = "complete.obs", nboot = 100)
# plot(fit1)
# pvrect(fit1, alpha=.8)

# x = seplot(fit, identify=FALSE)
# msplot(fit, edges=x)