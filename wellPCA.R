mainDir = "H:/Engineering_and_Hydro_Science/Projects/Groundwater_Resource_Assessment/groundwater_levels_NAVD88_Florida_Georgia/groundwater/"
new.data = "Y"

if(new.data == "Y"){
### FILLED ###
filled.Data = "Filled" # "Filled", Original" To use the filled data instead of the original data
if(filled.Data =="Original")(
  model.domain.data1 = read.csv(paste("C:/all.agency.mdata.csv", sep = ""), header = TRUE) ### all agency mdata
) else (model.domain.data1 = read.csv(paste(mainDir, "data/Filled_stations_stack_All_agency_MaxCorr_0.9_MatchPair_10_MinRegPeriod3.csv", sep = ""), header = TRUE)) ### all agency mdata
}
agency.mdata = model.domain.data1

clus_num = 3
cluster.subset = subset(GIS.inventory.cluster, GIS.inventory.cluster$cluster_num  == clus_num)
bad.station = as.character(cluster.subset$univName)
bad.station

## plot the station that shows as an anomalous cluster
stationName = bad.station
ind = match(agency.mdata$univName, stationName)
station.data = agency.mdata[!is.na(ind),]

start.date <- as.Date("01/01/2000",format = "%m/%d/%Y")  # Begin date for the analysis
end.date <- as.Date("12/31/2010",format = "%m/%d/%Y")  # End data for the analysis
station.data$monthYear = as.Date(station.data$monthYear, format = "%Y-%m-%d") # create data for monthYear
station.data.wide = dcast(station.data, monthYear~univName, value.var = "levelNAVD88")
station.data.wide = station.data.wide[which((station.data.wide$monthYear >= start.date) & 
                              (station.data.wide$monthYear <= end.date)),] # period of record 
st = na.omit(station.data.wide[,1:ncol(station.data.wide)])
sc.st = scale(st[,2:ncol(station.data.wide)])
pr = prcomp(sc.st)
screeplot(pr)

PC1 = scale(pr$x[,1])
plot(st$monthYear,PC1)
points(st$monthYear, sc.st[,1], col = 2)
points(st$monthYear, sc.st[,2], col = 3)
points(st$monthYear, sc.st[,3], col = 4)
points(st$monthYear, sc.st[,4], col = 5)
