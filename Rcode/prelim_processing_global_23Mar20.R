# processing of the climate data (already pre-processed by Rachel)
# see https://rfrelat.github.io/Spatial2_MultiExamples.html for nice script
# files on GitHub
# owner: Kath O'Reilly (LSHTM)
# Steps: bring in the data from Copernicus, WorldPop, and WHO
# for each country;
# 1) find temp & humidity over 01Jan-14Mar, adjusting for population
# 2) identify current WHO transmission status
# 3) couple of plots

setwd("~/Documents/GitHub/COVID19-season-short/Data")

#install.packages(c("maps", "mapdata", "ncdf4", "raster", "rgdal", "RColorBrewer", "sp"))

# To load and process GIS data
require(sp)
require(rgdal)
require(raster)
require(ncdf4)
#To make nicer looking maps
require(maps) 
require(mapdata)
require(RColorBrewer)
require(ptinpoly)
require(tmap)
require(sf)
require(ggplot2)
require(ggrepel)
library(grid)
library(snowfall)

# 0 - Load in data

dir <- "~/Documents/GitHub/COVID19-season-short/Shapefiles/Detailed/detailed_2013.shp" #"~/Dropbox/Covid-19-season/Rcode/Admin1(2011)/admin1.shp"
name <- "detailed_2013"
admn0 <- readOGR(dir,name)
shpdat0 <- data.frame(admn0@data)
dim(shpdat0)  # 224 polygons 237 for the detailed one
table(shpdat0$CNTRY_TERR)

# improve how the legend is plotted
#jpeg("test01.jpeg",height=1000,width=2000)
#tm_shape(admn0) + tm_polygons("WHO_REGION") + 
#  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
#            legend.title.size=2.2)
#dev.off()
#write.csv(shpdat0,file="~/Dropbox/Covid-19-season/Data/shpdat0B_20Mar20.csv",fileEncoding="UTF-8")  # UT8

# load in temperature files *updated*
tmean <- brick("~/Documents/GitHub/COVID19-season-short/Data/era5_daily_t2m_mean_brick_2019-11-01_2020-03-14.tif")
tmp <- tmean[[2]]
names(tmean)
dim(tmean)

#jpeg("test01.jpeg",height=2000,width=4000)
#plot(tmp)
#dev.off()

humid <- brick("~/Documents/GitHub/COVID19-season-short/Data/era5_relhum_mean_brick_2019-11-01_2020-03-14.tif")
tmp <- humid[[1]]
names(humid)
dim(humid)

#jpeg("test01.jpeg",height=2000,width=4000)
#plot(tmp)
#dev.off()

# load in pop raster
# WORLDPOP see https://www.worldpop.org/
str_name<-'~/Documents/GitHub/COVID19-season-short/Data/ppp_2020_1km_Aggregated.tif' 
popdat<-raster(str_name)
dim(popdat)

# now we want to 
# 1. select the relevant shp for each location
# 2. extract the relevant temp data 
# 3. extract the relevant pop data
# working with SHP areas to extract data for specific countries.

names(shpdat0)
shpdat0$temp_mean_14Mar20 <- shpdat0$temp_mean_14Mar20_pa <- 0  # temp on the 14Mar
shpdat0$temp_mean_3mth <- shpdat0$temp_mean_3mth_pa <- 0        # temp Jan-14Mar
shpdat0$humd_mean_14Mar20 <- shpdat0$humd_mean_14Mar20_pa <- 0  # temp on the 14Mar
shpdat0$humd_mean_3mth <- shpdat0$humd_mean_3mth_pa <- 0        # temp Jan-14Mar

dates <- seq(as.Date("2019/11/01"), as.Date("2020/03/14"), "days")  # use this sequence
# check we've got the right indexes for dates
dates[62] 
dates[135]
length(dates[122]:dates[31]) # taking summary of 92 days in total.

# loop through - this takes a while so next time should really have this in parallel
getRas_funct <- function(i,admn0,tmean,popdat,humid){
  tmp_shp <- admn0[i,]   # shp for country i
  # 1. raster of 14Mar20
  tmean01 <- tmean[[135]]
  craster01 <- crop(tmean01,tmp_shp)
  popras <- crop(popdat,tmp_shp)
  popras2 <- resample(popras,tmean01)  # scale to craster
  t1 <- as.vector(unlist(extract(popras2,tmp_shp))) # '1000
  t1[is.na(t1)] <- 0  # remove na's
  #tmean01
  # one temp for 14 March
  pop <- t1/sum(t1)
  temp <- as.vector(unlist(extract(tmean01,tmp_shp)))
  temp_mean_14Mar20 <- mean(temp) #cellStats(craster, stat='mean', na.rm=TRUE, asSample=TRUE)  # not pop adjusted
  temp_mean_14Mar20_pa <- sum(pop*temp)
  
  # humidity
  hmean01 <- humid[[135]]
  hraster01 <- crop(hmean01,tmp_shp)
  hvec <- as.vector(unlist(extract(hmean01,tmp_shp)))
  hmid_mean_14Mar20 <- mean(hvec) #cellStats(craster, stat='mean', na.rm=TRUE, asSample=TRUE)  # not pop adjusted
  hmid_mean_14Mar20_pa <- sum(pop*hvec)
  #browser()
  # mean temp for Dec-Feb
  temp_all <- c()
  humd_all <- c()
  for(t in 62:135){  # 2020-01-01" to "2020-03-14"
    # temp
    tmean01 <- tmean[[t]]   # extract the right raster
    temp <- as.vector(unlist(extract(tmean01,tmp_shp)))
    temp_all <- rbind(temp_all,temp)  # increasing rows
    # humidity
    hmean01 <- humid[[t]]   # extract the right raster
    hvec <- as.vector(unlist(extract(hmean01,tmp_shp)))
    humd_all <- rbind(humd_all,hvec)  # increasing rows
  }
  temp_mean_3mth <- mean(apply(temp_all,2,mean))
  temp_mean_3mth_pa <- sum(apply(temp_all,2,mean)*pop)       # temp 01Jan-14Mar
  hmid_mean_3mth <- mean(apply(humd_all,2,mean))
  hmid_mean_3mth_pa <- sum(apply(humd_all,2,mean)*pop)       # temp 01Jan-14Mar
  write.csv(paste0(as.character(shpdat0$CNTRY_TERR[i]," ",i," complete")),file=paste0("Checks/check_",i,".csv"))
  return(list=c(temp_mean_14Mar20=temp_mean_14Mar20,temp_mean_14Mar20_pa=temp_mean_14Mar20_pa,
                hmid_mean_14Mar20=hmid_mean_14Mar20,hmid_mean_14Mar20_pa=hmid_mean_14Mar20_pa,
                temp_mean_3mth=temp_mean_3mth,temp_mean_3mth_pa=temp_mean_3mth_pa,
                hmid_mean_3mth=hmid_mean_3mth,hmid_mean_3mth_pa=hmid_mean_3mth_pa))
}

#shpdat0[201,]
#tmp <- getRas_funct(201,admn0,tmean,popdat,humid)

cpuStart <- 1
cpuStop <- dim(shpdat0)[1] #total nodes
# considering above we need each node to specify a model
cpus=4
process <- function(parallel = FALSE){
  sfInit(parallel = parallel, cpus = cpus)
  # bespoke functions need to be exported
  # libraries need to be exported
  sfLibrary(raster)
  #sfLibrary(deSolve)
  sfExportAll()
  out <- sfLapply( cpuStart:cpuStop, getRas_funct,admn0=admn0,tmean=tmean,popdat=popdat,humid=humid)
  sfStop()
  return(out)
}

# ONLY RUN THIS IF YOU WANT TO WAIT 90MINS!
#in serial
#system.time(outS <- process(parallel = FALSE))
# in parallel
system.time(outP <- process(parallel = TRUE))
# save this
save(outP,file="temp_humd_out_detail_Jan-14Mar.Rdat")

# or load it here
load(file="~/Documents/GitHub/COVID19-season-short/Data/temp_humd_out_detail_Jan-14Mar.Rdat")
# library(data.table)
out <- data.frame(matrix(unlist(outP,use.names = TRUE),nrow=length(outP), byrow=T))
names(out) <- names(outP[[1]])

shpdat0$temp_mean_14Mar20 <- out$temp_mean_14Mar20
shpdat0$temp_mean_14Mar20_pa <- out$temp_mean_14Mar20_pa # temp on the 14Mar
shpdat0$temp_mean_3mth <- out$temp_mean_3mth
shpdat0$temp_mean_3mth_pa <- out$temp_mean_3mth_pa        # temp Dec/Jan/Feb
shpdat0$humd_mean_14Mar20 <- out[,3]
shpdat0$humd_mean_14Mar20_pa <- out[,4]  # temp on the 14Mar
shpdat0$humd_mean_3mth <- out[,7]
shpdat0$humd_mean_3mth_pa <- out[,8]       # temp Dec/Jan/Feb

#save(shpdat0,file="shpdat0_inctemphumid_updated_25Mar20.Rdat")

admn0$temp_mean_14Mar20 <- shpdat0$temp_mean_14Mar20
admn0$temp_mean_14Mar20_pa <- shpdat0$temp_mean_14Mar20_pa # temp on the 14Mar
admn0$temp_mean_3mth <- shpdat0$temp_mean_3mth
admn0$temp_mean_3mth_pa <- shpdat0$temp_mean_3mth_pa        # temp Dec/Jan/Feb
admn0$humd_mean_14Mar20 <- shpdat0$humd_mean_14Mar20
admn0$humd_mean_14Mar20_pa <- shpdat0$humd_mean_14Mar20_pa  # temp on the 14Mar
admn0$humd_mean_3mth <- shpdat0$humd_mean_3mth
admn0$humd_mean_3mth_pa <- shpdat0$humd_mean_3mth 

setwd("~/Documents/GitHub/COVID19-season-short/Plots")
# plot these to check we're outputting something sensible
jpeg("world_temp_14MarB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("temp_mean_14Mar20",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_temp_14Mar_paB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("temp_mean_14Mar20_pa",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_temp_3mthB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("temp_mean_3mth",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_temp_3mth_paB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("temp_mean_3mth_pa",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_humd_14MarB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("humd_mean_14Mar20",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_humd_14Mar_paB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("humd_mean_14Mar20_pa",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_humd_3mthB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("humd_mean_3mth",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()
jpeg("world_humd_3mth_paB.jpeg",height=1000,width=2000)
tm_shape(admn0) + tm_polygons("humd_mean_3mth_pa",palette="PRGn") + 
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=2.2)
dev.off()

# Add WHO data

sitrep <- read.csv("~/Documents/GitHub/COVID19-season-short/Data/WHO_SITREP_COVID_24032020_v2.csv")
names(sitrep)
# ISO3 has been added in by hand to facilitate matching

shpdat0$cases <-  shpdat0$new_cases <- shpdat0$deaths <- shpdat0$new_deaths <- shpdat0$transmission <- shpdat0$days_since_last_case <- "No data" 
oo <- match(as.character(shpdat0$ISO_3_CODE),as.character(sitrep$ISO_3_CODE))
sum(!is.na(oo))
sum(is.na(oo))

shpdat0$cases <- sitrep$Total.confirmed[oo]
shpdat0$new_cases <- sitrep$Total.confirmed.new.cases[oo]
shpdat0$deaths <- sitrep$Total.deaths[oo]
shpdat0$new_deaths <- sitrep$Total.new.deaths[oo]
shpdat0$transmission <- sitrep$Transmission.classification[oo]
shpdat0$days_since_last_case <- sitrep$Days.since.last.reported.case[oo]

# cleaning
# ember dat from... the Diamond princess data is included in sitreps - ignore.
table(sitrep$Total.deaths)
table(shpdat0$deaths) 

table(sitrep$Transmission.classification)
table(shpdat0$transmission)
shpdat0$transmission <- relevel(shpdat0$transmission,levels(shpdat0$transmission)[2])
levels(shpdat0$transmission)

# into maps
admn0$transmission <- shpdat0$transmission
levels(admn0$transmission) <- c("Local transmission ","Imported cases only ","Under investigation ","No information")
admn0$transmission[is.na(admn0$transmission)] <- "No information"
admn0$deathsGrp <- shpdat0$deathsGrp
table(admn0$transmission)

#####################################################
# Mapping to shpaefiles countries
#####################################################

# choose a colour-scheme - using https://colorbrewer2.org/ 
myColors <- c("#998ec3","#f1a340","#f5f5f5","grey80")
sc <- 2

jpeg("sitrep_WHOdef_24Mar20v2.jpeg",height=450*sc,width=900*sc)
tm_shape(admn0) + tm_polygons(title="WHO transmission status\n as of 24 Mar 2020","transmission",palette=myColors) +
  tm_layout(legend.position=c("left","bottom"),legend.text.size=2,  # ok this is being recognised
            legend.title.size=4)
dev.off()

# so now... I think it's better to link directly to sitrep rather than other way round 
# (the other way round is too messy)
# dim(sitrep)

head(sitrep)
oo <- match(as.character(sitrep$ISO_3_CODE),as.character(shpdat0$ISO_3_CODE))
sum(is.na(oo))

sitrep$temp_mean_14Mar20 <- shpdat0$temp_mean_14Mar20[oo]
sitrep$temp_mean_14Mar20_pa <- shpdat0$temp_mean_14Mar20_pa[oo] # temp on the 14Mar - populatin adjusted
sitrep$temp_mean_3mth <- shpdat0$temp_mean_3mth[oo]
sitrep$temp_mean_3mth_pa <- shpdat0$temp_mean_3mth_pa[oo]        # temp 01Jan-14Mar
sitrep$humd_mean_14Mar20 <- shpdat0$humd_mean_14Mar20[oo]
sitrep$humd_mean_14Mar20_pa <- shpdat0$humd_mean_14Mar20_pa[oo]  # temp on the 14Mar
sitrep$humd_mean_3mth <- shpdat0$humd_mean_3mth[oo]              # 01Jan-14Mar
sitrep$humd_mean_3mth_pa <- shpdat0$humd_mean_3mth[oo] 

sitrep[is.na(oo),]
# ok so we don't have infomration for Curcao, Guensey, Jersey, Saint Barth, Saint Martin.
# I'm justifying the exclusion of these because they are small islands (and/or part of UK)

# ok so each of these need to be identified...
sitreptmp <- sitrep[is.na(oo),]
sitreptmp <- sitrep[is.na(sitrep$temp_mean_14Mar20_pa),]  # 19 obs - all islands (Singapore being the largest)

# make small adjustment for pa being NA (small number of countries)
sitrep$temp_mean_14Mar20_pa[is.na(sitrep$temp_mean_14Mar20_pa)] <- sitrep$temp_mean_14Mar20[is.na(sitrep$temp_mean_14Mar20_pa)]
sitrep$temp_mean_14Mar20_paadj <- sitrep$temp_mean_14Mar20_pa
sitrep$temp_mean_3mth_pa[is.na(sitrep$temp_mean_3mth_pa)] <- sitrep$temp_mean_3mth[is.na(sitrep$temp_mean_3mth_pa)]
sitrep$temp_mean_3mth_paadj <- sitrep$temp_mean_3mth_pa

# fileEncoding hopefully gets ride of the word gremlins
write.csv(sitrep,file="~/Documents/GitHub/COVID19-season-short/Data/WHO_SITREP_COVID_24032020_incTemp_v2.csv",fileEncoding="UTF-8")

# so let's focus on 'local transmission' and 'importations'
sitrepB <- sitrep[sitrep$Transmission.classification!="Under investigation",]

labdea <- sitrep[nchar(sitrep$labs)>1,]

sitrepB$labs <- " "
# sitrepB$labs[sitrepB$Transmission.classification==levels(sitrepB$Transmission.classification)[2] & 
#                sitrepB$temp_mean_3mth_pa>25 & 
#                sitrepB$humd_mean_3mth_pa>70 &
#                !is.na(sitrepB$temp_mean_3mth_pa) & 
#                !is.na(sitrepB$Transmission.classification)] <- 
#   as.character(sitrepB$Reporting.Country[sitrepB$Transmission.classification==levels(sitrepB$Transmission.classification)[2] & 
#                                            sitrepB$temp_mean_3mth_pa>25 & 
#                                            sitrepB$humd_mean_3mth_pa>70 &
#                                            !is.na(sitrepB$temp_mean_3mth_pa) & 
#                                            !is.na(sitrepB$Transmission.classification)])
sitrepB$labs[sitrepB$Transmission.classification==levels(sitrepB$Transmission.classification)[2] & 
               sitrepB$temp_mean_3mth_pa>25 & 
               sitrepB$humd_mean_3mth_pa>70 &
               !is.na(sitrepB$temp_mean_3mth_pa) & 
               !is.na(sitrepB$Transmission.classification)] <- "                           +"

# make dataframe just for labels
labtrs <- sitrepB[nchar(as.character(sitrepB$labs))>1,]
table(as.character(labtrs$Reporting.Country))

myColors <- c("#f1a340","#998ec3","#f5f5f5","grey80")
jpeg("compare_temp3mth_popad_WHOdef_24Mar20v2c.jpeg",height=400*sc,width=900*sc*0.5)
p1 <- ggplot(sitrepB,aes(x=Transmission.classification,y=temp_mean_3mth_pa,
                         colours=Transmission.classification)) + #geom_point() +
  geom_jitter(aes(colour = factor(Transmission.classification)),shape=16,size=6,
              position=position_jitter(0.2),show.legend = FALSE) +
  scale_colour_manual(name="test",values = myColors) + #geom_text(data=sitrep,label=labs) + 
  labs(y="Mean temp (oC pop adj.)",x="WHO definition (as of 24 Mar 2020)",size=6) +
  theme_bw(base_size = 40)
p1 + geom_text(data=labtrs,aes(label=labs),#nudge_x = 0.5,#check_overlap = T,#direction="both",
               position=position_jitter(width=0.1,height=0),
                     segment.colour = NA,size=12)
dev.off()
# notes:
# Mogolia is -15 and local transmission only.

# and humidity
jpeg("compare_humid3mth_popad_WHOdef_24Mar20v2b.jpeg",height=400*sc,width=900*sc*0.5)
p1 <- ggplot(sitrepB,aes(x=Transmission.classification,y=humd_mean_3mth_pa,
                         colours=Transmission.classification)) + #geom_point() +
  geom_jitter(aes(colour = factor(Transmission.classification)),shape=16,size=6,
              position=position_jitter(0.2),show.legend = FALSE) +
  scale_colour_manual(name="test",values = myColors) + #geom_text(data=sitrep,label=labs) + 
  labs(y="Mean rel humidity (% pop adj.)",x="WHO definition (as of 24 Mar 2020)",size=6) +
  theme_bw(base_size = 40)
p1 + geom_text(data=labtrs,aes(label=labs),#nudge_x = 0.5,#check_overlap = T,#direction="both",
               position=position_jitter(width=0.1,height=0),
               segment.colour = NA,size=12)
dev.off()
# notes:
# BFA (Burkina Faso) has <20% humidity and local transmission



# end