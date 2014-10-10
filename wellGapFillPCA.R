#########
# TITLE: PRINCIPAL COMPONENT FILLING FOR GROUNDWATER WELLS 
# AUTHOR: NATHAN JOHNSON
# DATE: 2014-06-27

#packages for data formatting
require(reshape2)
require(stats)
require(gdata)
require(graphics)
require(grDevices)
require(utils)
require(lubridate) # year() function
require(maptools)
require(maps)
require(mapdata)
require(mapproj)
require(plyr) # count() function

#### USER INPUT ####
# station.list = "USGSGA314513083104101"

# pc = "PC_14"
# max.lat = 32.25
# min.lat = 31
# min.lon = -82.5
# max.lon = -81.5
# station.exclude = c("USGSGA320452082201401","USGSGA311322082052801","USGSGA321239082060301", "USGSGA310706082155101", "USGSGA320436082185801", "USGSGA321100081492701", "USGSGA320452082071001", "USGSGA321032081535001","USGSGA320050082241501",  "USGSGA310320082161801", "USGSGA321219082050301", "USGSGA321415082174301", "USGSGA321240082231101", "USGSGA313936081502501")

# west
# pc = "PC_12"
# max.lat = 32.25
# min.lat = 31
# min.lon = -84
# max.lon = -82
# station.exclude = c("USGSGA311322082052801", "USGSGA321239082060301", "USGSGA321302082243601", "USGSGA321415082174301", "USGSGA321240082231101", "USGSGA322530081465601")

# north 
pc = "PC_12"
max.lat = 33
min.lat = 32
min.lon = -83
max.lon = -80.8
station.exclude = c("USGSGA311322082052801", "USGSGA321239082060301", "USGSGA321302082243601", "USGSGA321415082174301", "USGSGA321240082231101", "USGSGA322530081465601")

# pc = "PC_1"
# max.lat = 31.5
# min.lat = 30.25
# min.lon = -83
# max.lon = -82
# station.exclude = c("USGSGA310620082342201", "USGSGA311322082052801")

# pc = "PC_20"
# max.lat = 32.5
# min.lat = 32
# min.lon = -81.5
# max.lon = -80.5
# station.exclude = c("")

obs.min =  4 # minimum data available over the period selected
obs.max = 400 # 364 # maximum data available over the period selected
start.date <- as.Date("01/01/1982",format = "%m/%d/%Y")  # Begin date for the analysis
end.date <- as.Date("12/01/2012",format = "%m/%d/%Y")  # End data for the analysis

#### DATA DIRECTORY ####
mainDir = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/"
model.name = "GA" # which region for cluster analysis (CSM, NFSEG, All_agency)
aquifer.type = c("UFA") # c("FAS", "SAS")
filled.Data = "Filled" # "Filled", "Original"
new.data = "Y"

## Data
GIS.inventory = read.csv(paste(mainDir, "data/inventory/All_agency_well_inventory_strat.csv", sep =""))
if(new.data == "Y"){
  if(filled.Data =="Original"){
    model.data = read.csv(paste("C:/all.agency.mdata.csv", sep = ""), header = TRUE) ### all agency mdata
  } else {model.data = read.csv(paste(mainDir, "data/Filled_stations_stack_second_All_agency_MaxCorr_0.9_MatchPair_10_MinRegPeriod3.csv", sep = ""), header = TRUE)} ### all agency mdata
} else {model.data = model.domain.data[2:5]}

## Reduce data to only include model domain ##
if(model.name != "All_agency"){
  subset.inventory = subset(GIS.inventory, !is.na(GIS.inventory[model.name]))
  ind = match(model.data$univName, GIS.inventory$univName)
  subset.data = model.data[!is.na(ind),] 
} else {subset.data = model.data}

## Reduce data to exclude stations that are duplicates of other 
if(length(station.exclude) != 0){ 
  subset.inventory = subset(subset.inventory, !subset.inventory$univName %in% station.exclude)
  subset.data = subset(subset.data, subset.data$univName %in% subset.inventory$univName)
}

## Reduce data to only include aquifer.type ##
if(aquifer.type != ""){ 
  subset.inventory = subset(subset.inventory, subset.inventory$aquiferFinal %in% aquifer.type)
  subset.data = subset(subset.data, subset.data$univName %in% subset.inventory$univName)
}

## Reduce data by start and end dates ##
subset.data = subset(subset.data, subset = as.Date(subset.data$monthYear) >= start.date & as.Date(subset.data$monthYear) <= end.date)

## Reduce data by number of observations ##
data.avail = count(subset.data, "univName")
subset.inventory = merge(x=subset.inventory, y = data.avail, by = "univName")
subset.inventory = subset(subset.inventory, subset = subset.inventory$freq >= obs.min & 
                            subset.inventory$freq <= obs.max)
subset.data = merge(x=subset.inventory, y=subset.data, by = "univName")
subset.data = subset(subset.data, select = c("univName","monthYear", "levelNAVD88", "data.type"))

## Reduce data to spatial extent ##
subset.inventory =  subset(subset.inventory, min.lat < subset.inventory$latDD
                           & max.lat > subset.inventory$latDD
                           & min.lon < subset.inventory$lonDD
                           & max.lon > subset.inventory$lonDD)
subset.data = merge(x=subset.inventory, y=subset.data, by = "univName")
subset.data = subset(subset.data, select = c("univName","monthYear", "levelNAVD88", "data.type"))

####  DATA FORMATTING AND PERIOD SELECTION ####
wide.data = dcast(subset.data, monthYear~univName, value.var = "levelNAVD88")
wide.data$monthYear = as.Date(wide.data$monthYear, format = "%Y-%m-%d") # create data for monthYear
row.names(wide.data) = wide.data$monthYear
wide.data = wide.data[order(wide.data$monthYear),] # order the dataframe
wide.data = wide.data[which((wide.data$monthYear >= start.date) & 
                              (wide.data$monthYear <= end.date)),]

# Create several variables for the nested loops that follow #
stat.model.merge = data.frame(monthYear = character(), level = numeric(), stringsAsFactors = FALSE)
# stat.model.merge.all.wide = data.frame(monthYear = data.frame.reg$monthYear)
stat.model.merge.all.stack = data.frame(monthYear = as.Date(character()), 
                                        levelNAVD88 = numeric(), univName = character(), 
                                        data.type = character())
sum.stat.table = data.frame(orig_station = character(), 
                            fill_station = character(), 
                            R2 = numeric(), 
                            DF = numeric(), 
                            RMSE= numeric(),
                            intercept = numeric(),
                            slope = numeric(),
                            slope.SE = numeric())

station.list = as.character(unique(subset.inventory$univName))

for(st in 2:length(station.list)){
  station = colnames(wide.data[st])
  station.inventory = subset(subset.inventory, subset = subset.inventory$univName == station)
  station.fill = wide.data[c("monthYear", station)] #create small dataframe with only station and pc
  # prin.comp.fill = wide.data[c("monthYear", pc)]
  prin.comp.fill = prin.comp.all[c("date", pc)]
  names(prin.comp.fill)[1] = "monthYear"
  data.frame.reg = merge(station.fill, prin.comp.fill, by = "monthYear", all = TRUE)
  matched.points = nobs(data.frame.reg[[station]] + data.frame.reg[[pc]])
  if(matched.points <=1)(next())
  reg = lm(data.frame.reg[[station]]~data.frame.reg[[pc]])
  
  rsquare = data.frame(r.squared = numeric(0), df = numeric(0), 
                       station = character(0), stringsAsFactors = FALSE)  #create data.frame for R2 values
  sum = summary(reg)
  rsquare = rbind(rsquare, data.frame(r.squared = sum$r.squared, # dataframe of R2,DF for all stations
                                      df = sum$df[2],station = station))
  if(rsquare$r.squared < 0.7|reg$coefficients[2] < 0){next}
  #   stat.model.merge.all.stack = rbind(
  #     stat.model.merge.all.stack, 
  #     data.frame(monthYear = as.character(data.frame.reg$monthYear),
  #                levelNAVD88 = data.frame.reg[[station]], 
  #                univName = station, 
  #                data.type = c("original"), stringsAsFactors = FALSE))
  
  stat.model = predict.lm(object = reg, 
                          newdata = data.frame(data.frame.reg[[station]])) # create lm statistical model of subset.label using the best.station
  station.fill = data.frame(stat.model, data.frame.reg[[station]]) # filled station and stat.model dataframe
  ind = station.fill[,2]+station.fill[,1] # index of points that stat model will fill of original station
  ind1 = data.frame(data.type = as.character(ind), filled ="PCfilled", 
                    original = "original", stringsAsFactors = FALSE)
  ind1[is.na(ind),1] = ind1[is.na(ind),2]
  ind1[!is.na(ind),1] = ind1[!is.na(ind),3]
  station.fill[!is.na(ind),2] = station.fill[!is.na(ind),1] # create a dataframe with filled station and stat model
  
  stat.model.merge = data.frame(data.frame.reg$monthYear, station.fill[,1]) # station filled in with best regressed station.
  colnames(stat.model.merge) = c("monthYear", paste(station, "_filled_station", sep = ""))
  #   stat.model.merge.all.wide = cbind(stat.model.merge.all.wide, stat.model.merge[2])
  colnames(stat.model.merge) = c("monthYear", "levelNAVD88")
  
  stat.model.merge.all.stack = rbind(
    stat.model.merge.all.stack, 
    data.frame(monthYear = as.character(stat.model.merge$monthYear), 
               levelNAVD88 = stat.model.merge$levelNAVD88, 
               univName = station, data.type = ind1[1]))
  stat.model.merge.all.stack = na.omit(stat.model.merge.all.stack)
  sum.stat = data.frame(orig_station = station, 
                        fill_station = pc, 
                        R2 = format(summary(reg)$r.squared, digits = 3), 
                        DF = format(summary(reg)$df[2]), 
                        RMSE= format(sqrt(mean(reg$residuals^2)), digits = 3), 
                        intercept = reg$coefficients[1], 
                        slope = reg$coefficients[2], 
                        slope.SE = summary(reg)$coefficients[2,2])
  sum.stat.table = rbind(sum.stat.table, sum.stat)
  
  ##### FIGURES #####
  png(file = paste(mainDir, "figures/PCAgapFill/", station,"_", pc, "_PCAgapFill.png", sep = ""),
      width = 1000, height = 1000, units= "px", res=100)
  layout(matrix(c(5,5,3,3,5,5,3,3,1,2,4,4,1,2,4,4),4,4))
  # layout.show(nf)
  
  ## Hydrographs of station
  plot(data.frame.reg$monthYear,data.frame.reg[[station]], type = "p",
       main= paste("Original station: ", station, sep = ""), 
       pch = 20, xlab = "Date", ylab = "Level (NAVD88,ft)")
  
  ## Hydrograph of Principal Component
  plot(data.frame.reg$monthYear,data.frame.reg[[pc]], type = "p",
       main= paste("Principal Component: ", pc, sep = ""), 
       pch = 20, xlab = "Date", ylab = "Stnd_")
  
  ## Correlation firgure
  plot(data.frame.reg[[station]]~data.frame.reg[[pc]], 
       xlab = pc, ylab = station, main = "Regression")
  abline(reg)
  legend(x = "topleft", bty= "n",
         legend = c(paste("R2 =", format(summary(reg)$r.squared, digits = 3)), 
                    paste("RMSE =", format(sqrt(mean(reg$residuals^2)), digits = 3)),
                    paste("DF =", format(summary(reg)$df[2]))))
  
  ## Filled station
  plot(data.frame.reg$monthYear[is.na(ind)], 
       station.fill[is.na(ind),1], col = "red", type = "p", pch = 1, 
       main = paste("Filled station:", station), ylab = "Level (NAVD88,ft)", xlab = "Date", 
       xlim = c(start.date, end.date), 
       ylim = c(min(station.fill[is.na(ind),1],data.frame.reg[[station]], na.rm= TRUE), max(station.fill[is.na(ind),1],data.frame.reg[[station]], na.rm = TRUE))
  )
  points(data.frame.reg$monthYear, data.frame.reg[[station]],
         col = "black", type = "p", pch = 20, cex = 2)
  legend("bottomleft", bty = "n", c("Original data", "PC Filled data"), col = c("black","red"), 
         lty = c(0, 0),pch=c(20,1), merge = F, border=F, cex = 1)
  
  # Map of well location
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
  points(x = station.inventory$lonDD, y = station.inventory$latDD, col = 'red', pch = 3, cex = 1)
  text(x = station.inventory$lonDD, y = station.inventory$latDD, 
       labels = station.inventory$univName, col = "black", pos = 4, cex = 1)
  dev.off()
}
PCfilled.data = subset(stat.model.merge.all.stack, subset = stat.model.merge.all.stack$data.type == "PCfilled")

# write.table(PCfilled.data, file = paste(mainDir, "data/PCAgapfill/PCfilled.csv", sep = ""), append = TRUE, row.names = FALSE,col.names = FALSE, sep = ",")
# write.table(sum.stat.table, file = paste(mainDir, "data/PCAgapfill/PCfilled_Stat_table.csv", sep = ""), append = TRUE, row.names = FALSE, sep = ",",col.names = FALSE)
# 
# # ## CREATE DATA TABLES FOR GIS WITH ALL THE MONTHLY AND ANNUAL VALUES
# all.data = model.domain.data1[,2:5]
# PC.data = read.csv(file = paste(mainDir, "data/PCAgapfill/PCfilled.csv",sep = ""), header = FALSE)
# colnames(PC.data) = names(PCfilled.data)
# all.data = rbind(all.data, PC.data)
# all.data$monthYear = as.Date(all.data$monthYear)
# all.data = data.frame(all.data, year = paste("WL_",year(all.data$monthYear), sep =""))
# all.data.avg = data.frame(all.data[1:3], data.type = "Avg", all.data[5])
# all.data = rbind(all.data, all.data.avg)
# concat.df = data.frame(all.data, concat = paste(all.data$univName, all.data$data.type, sep = "_"))
# concat.df$concat = sub(x = concat.df$concat, pattern = "filled_second", replacement = "filled.second", fixed = TRUE)
# 
# #monthly data tables
# wide.all.data.month = dcast(concat.df, concat~monthYear, value.var = "levelNAVD88", fun.aggregate = mean, fill = -9999) # create wide data.frame with concat.df
# names(wide.all.data.month)[2:length(wide.all.data.month)] = paste("WL_",names(wide.all.data.month[2:length(wide.all.data.month)]), sep = "")
# name.type = data.frame(do.call('rbind', strsplit(as.character(wide.all.data.month$concat), split = '_', fixed = TRUE)))
# wide.all.data.month.final = cbind(name.type, wide.all.data.month)
# names(wide.all.data.month.final)[1:2] = c("univName", "data.type")
# 
# #annual data tables
# wide.all.data.year = dcast(concat.df, concat~year, value.var = "levelNAVD88", fun.aggregate = median, fill = -9999)
# name.type = data.frame(do.call('rbind', strsplit(as.character(wide.all.data.year$concat), split = '_', fixed = TRUE)))
# wide.all.data.year.final = cbind(name.type, wide.all.data.year)
# names(wide.all.data.year.final)[1:2] = c("univName", "data.type")
# wide.all.data.year.final1 = wide.all.data.year.final
# wide.all.data.year.final1[,4:length(wide.all.data.year.final)] = as.numeric(wide.all.data.year.final[,4:length(wide.all.data.year.final)])
# apply
# 
# ## WRITE DATA TO .CSV FILE
# write.table(wide.all.data.month.final, file = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/data/GISlevelsMonth.csv", row.names =  FALSE, sep = ",")
# write.table(wide.all.data.year.final, file = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/data/GISlevelsYear.csv", row.names =  FALSE, sep = ",")

# library(xlsx)
# write.xlsx(wide.all.data.month.final, file = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/data/GISlevelsMonth.xlsx", row.names =  FALSE, sep = ",", append = FALSE, showNA = TRUE, col.names = TRUE, sheetName = "GIS_WL_month")
# write.xlsx(wide.all.data.year.final, file = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/data/GISlevelsYear.xlsx", row.names =  FALSE, sep = ",", append = FALSE, showNA = TRUE, col.names = TRUE, sheetName = "GIS_WL_month")