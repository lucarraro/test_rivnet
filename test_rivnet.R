# test_rivnet ####
# Scripts supporting Carraro, L., "Seamless extraction and analysis of river networks in R:  the rivnet package"
# NOTE: to produce Fig. 5, you need to download the landcover map of Switzerland (see lines 238-239)

rm(list=ls())
library(rivnet) # rivnet >=0.3.0 must be used
library(sp)
library(rnaturalearth)
library(terra)

# FIG 2 ####
rivers <- read.csv("rivers4.csv")
pdf("Fig2.pdf", width = 21/2.54, height = 12/2.54)
par(mfrow=c(1,2))
# extract river with reduced zoom to show plot (river lines appear thicker)
river <- rivnet::extract_river(outlet = c(rivers$x_outlet[1], rivers$y_outlet[1]),
                               EPSG = rivers$EPSG[1],
                               ext = c(rivers$x_ll[1], rivers$x_tr[1], rivers$y_ll[1], rivers$y_tr[1]),
                               z = rivers$zoom[1] - 1,
                               displayUpdates = 0,
                               showPlot = TRUE)
# extract same river with higher resolution to show functioning of locate_site
river <- rivnet::extract_river(outlet = c(rivers$x_outlet[1], rivers$y_outlet[1]),
                               EPSG = rivers$EPSG[1],
                               ext = c(rivers$x_ll[1], rivers$x_tr[1], rivers$y_ll[1], rivers$y_tr[1]),
                               z = rivers$zoom[1],
                               displayUpdates = 0)
river <- rivnet::aggregate_river(river, thrA = 1e6)
locate_site(410000, 5200000, river, showPlot = TRUE)
dev.off()

# FIGURE 3 ####
rivers <- read.csv("rivers4.csv")
if (!file.exists("results_4rivers.rda")){
  out_L <- list(length(rivers$Name))
  out_A <- list(length(rivers$Name))
  out_PA <- list(length(rivers$Name))
  out_dA <- list(length(rivers$Name))
  river_polygons <- list(length(rivers$Name))
  computationalTime <- matrix(0,5,4)
  cell_ratio <- no_cells <- numeric(4)
  for (i in 1:length(rivers$Name)){
    t0 <- Sys.time()
    river <- rivnet::extract_river(outlet = c(rivers$x_outlet[i], rivers$y_outlet[i]),
                                   EPSG = rivers$EPSG[i],
                                   ext = c(rivers$x_ll[i], rivers$x_tr[i], rivers$y_ll[i], rivers$y_tr[i]),
                                   z = rivers$zoom[i],
                                   displayUpdates = 1,
                                   n_processes = 8,
                                   showPlot = TRUE)
    no_cells[i] <- river$FD$nNodes
    cell_ratio[i] <- river$FD$nNodes/river$dimX/river$dimY
    t1 <- Sys.time()
    river <- rivnet::aggregate_river(river, thrA = 1e6, displayUpdates = 1)
    t2 <- Sys.time()
    river <- rivnet::paths_river(river, level = "AG", displayUpdates = 1)
    t3 <- Sys.time()

    river_polygons[[i]] <- Polygon(cbind(river$CM$XContour[[1]][[1]], river$CM$YContour[[1]][[1]]))

    # max length upstream
    channelHeads <- which(river$AG$nUpstream==1) # find headwaters
    L <- Lhead <- numeric(river$AG$nNodes)
    for (j in channelHeads){ Lhead[j] <- max(sqrt((river$FD$X[river$SC$toFD[[j]]] - river$AG$X[j])^2 +
                                                    (river$FD$Y[river$SC$toFD[[j]]] - river$AG$Y[j])^2))}

    for (j in 1:river$AG$nNodes){ # add Lhead only to sources upstream of j
      if (!(j %in% channelHeads)){L[j] <- max(river$AG$downstreamPathLength[,j] + Lhead*(river$AG$downstreamPathLength[,j]>0))}}

    L <- L[-channelHeads] # only pick confluences
    A <- river$AG$A[-channelHeads]

    out_L[[i]] <- L
    out_A[[i]] <- A

    t4 <- Sys.time()

    # probability of aggregation
    dA <- sort(unique(unique(river$FD$A)))
    PA <- numeric(length(dA))
    for (j in 1:length(dA)){PA[j] <- sum(river$FD$A>=dA[j])/river$FD$nNodes }

    out_PA[[i]] <- PA
    out_dA[[i]] <- dA
    t5 <- Sys.time()
    computationalTime[1, i] <- difftime(t1,t0,unit="secs")
    computationalTime[2, i] <- difftime(t2,t1,unit="secs")
    computationalTime[3, i] <- difftime(t3,t2,unit="secs")
    computationalTime[4, i] <- difftime(t4,t3,unit="secs")
    computationalTime[5, i] <- difftime(t5,t4,unit="secs")
    save(out_L,out_A,out_PA,out_dA,river_polygons,computationalTime,
         no_cells,cell_ratio,
         file="results_4rivers.rda")

  }} else {(load("results_4rivers.rda"))}


kol <- hcl.colors(4,"PiYG")

pdf("Fig3.pdf", width = 21/2.54, height = 12/2.54)
h_values <- numeric(4)
par(mfrow = c(1,2))
# Hack's law plot
for (i in seq(4,1,-1)){
  L <- out_L[[i]]; A <- out_A[[i]]
  if (i==4){
    plot(A/1e6, L/1e3, col = kol[i], pch = 19, cex = 0.4,log = "xy",
         xlab = "A [km2]", ylab = "L [km]", xlim = c(1, 1e4), ylim = c(1, 200))
  } else {points(A/1e6, L/1e3, col = kol[i], pch = 19, cex = 0.4)}
  lmod <- lm(log(L/1e3) ~ log(A/1e6))
  ss <- summary(lmod); h_values[i] <- ss$coefficients[2,1]
  xx <- pracma::logspace(log10(min(A/1e6)), log10(max(A/1e6)), 20)
  lines(xx, exp(ss$coefficients[1,1])*xx^ss$coefficients[2,1], col = kol[4], lwd = 0.5)
}
legend(1, 200, legend = c(sprintf("Ilfis  -  h = % .3f",h_values[1]),
                          sprintf("Toss  -  h = % .3f", h_values[2]),
                          sprintf("Inn  -  h = % .3f", h_values[3]),
                          sprintf("Rhone  -  h = % .3f", h_values[4])),
       col = kol, pch = 19, text.font = 1, pt.cex = 0.5)

# scaling of drainage areas
cex_vec <- seq(0.2, 0.8, 0.2); thr <- 0.025
beta_values <- numeric(4)
for (i in seq(4,1,-1)){
  dA <- out_dA[[i]]; PA <- out_PA[[i]]
  if (i==4){
    plot(dA/1e6, PA, col = kol[i], pch = 19, cex = cex_vec[i], log = "xy",
         xlim = c(1e-4, 1e4), ylim = c(1e-8, 2), xlab = "a [km2]", ylab = "P[A >= a]")
  } else {points(dA/1e6, PA, col = kol[i], pch = 19, cex = cex_vec[i])}
  lmod <- lm(log(PA[dA <= thr*max(dA)]) ~ log(dA[dA <= thr*max(dA)]/1e6)) # cutoff: 5% of max A
  ss<- summary(lmod); beta_values[i] <- ss$coefficients[2,1]
  xx <- pracma::logspace(log10(min(dA/1e6)), log10(max(dA/1e6)), 20)
  lines(xx, exp(ss$coefficients[1,1])*xx^ss$coefficients[2,1], col = kol[4], lwd = 0.5)
}
legend(1e-3, 1e-3, legend = c(sprintf("Ilfis  -  beta = % .3f",-beta_values[1]),
                              sprintf("Toss  -  beta = % .3f", -beta_values[2]),
                              sprintf("Inn  -  beta = % .3f", -beta_values[3]),
                              sprintf("Rhone  -  beta = % .3f", -beta_values[4])),
       col = kol, pch = 19, text.font = 1, pt.cex = cex_vec)
dev.off()

# MAP FIGURE
pdf("Fig3_inset.pdf", width = 21/2.54, height = 12/2.54)
CH <- ne_countries(10, country = "Switzerland")
multicountry <- ne_countries(10, country = c("Switzerland",'Italy',"France","Germany","Austria"))
plot(multicountry, xlim=bbox(CH)[1,], ylim=bbox(CH)[2,])
abline(h=c(46,46.5,47,47.5), v=c(6,8,10),col="gray80")
for (i in 1:length(rivers$Name)){
  ps <- Polygons(list(river_polygons[[i]]), 1)
  xy <- SpatialPolygons(list(ps), proj4string = CRS(paste0("+init=epsg:23032")))
  XY <- spTransform(xy, CRS("+init=epsg:4326"))
  plot(XY, add = T,col = kol[i])
}
# add Thur
Thur <- rivnet::extract_river(outlet = c(517430, 5259917),
                              EPSG = 23032,
                              ext = c(481616, 552566, 5219109, 5267687),
                              z = 9)
Thur_polygon <- Polygon(cbind(Thur$CM$XContour[[1]][[1]], Thur$CM$YContour[[1]][[1]]))
ps <- Polygons(list(Thur_polygon),1)
xy <- SpatialPolygons(list(ps), proj4string = CRS(paste0("+init=epsg:23032")))
XY <- spTransform(xy, CRS("+init=epsg:4326"))
plot(XY, add = T,col = "grey50")
#abline(h = seq(46, 47.5, 0.5), v = 6:10, col = "grey90")

legend(x = 6, y = 48, legend = c(rivers$Name, "Thur"), fill = c(kol,"grey50"))
dev.off()

# FIGURE 4 ####
time_breakdown <- read.csv("times.csv")
time_elevatr <- matrix(as.numeric(time_breakdown[1,2:21]),ncol=4,nrow=5,byrow=TRUE)
time_tauDEM <- matrix(as.numeric(time_breakdown[2,2:21]),ncol=4,nrow=5,byrow=TRUE)
time_extract <- matrix(as.numeric(time_breakdown[3,2:21]),ncol=4,nrow=5,byrow=TRUE)
time_aggregate <- matrix(as.numeric(time_breakdown[4,2:21]),ncol=4,nrow=5,byrow=TRUE)
time_paths <- matrix(as.numeric(time_breakdown[5,2:21]),ncol=4,nrow=5,byrow=TRUE)

tot_time <- (colMeans(time_elevatr)+colMeans(time_tauDEM)+
               colMeans(time_extract)+colMeans(time_aggregate)+colMeans(time_paths))

breakdown <- rbind(colMeans(time_elevatr),colMeans(time_tauDEM),
                   colMeans(time_extract),colMeans(time_aggregate),colMeans(time_paths))
for (i in 1:length(tot_time)) breakdown[,i] <- breakdown[,i]/tot_time[i]

pdf("Fig4.pdf", width = 21/2.54, height = 10/2.54)
kols <- hcl.colors(5,"Zissou 1")
par(mfrow=c(1,2))
plot(no_cells,tot_time,log="x", ylim=c(0,60),xlim=c(5e4,2e6),
     bty="n",xaxt="n",yaxt="n",type="b",pch=19,ylab="Time [s]",xlab="Number of cells")
axis(1,pos=0); axis(2,pos=5e4)
barplot(breakdown,col=kols,names.arg=c("Ilfis","Toss","Inn","Rhone"))
dev.off()


# FIGURE 5 ####
Thur <- rivnet::extract_river(outlet = c(735010, 261530),
                              EPSG = 21781,
                              ext = c(700000, 770000, 220000, 270000),
                              z = 9)

Thur <- rivnet::aggregate_river(Thur, thrA = 1e6, maxReachLength = 2500)

# Hydrological stations data can be retreived from https://www.hydrodaten.admin.ch/en/stations-and-data.html
# legend: Thur (ID 2303) - Rietholzbach (ID 2414) - Necker (ID 2374) - Glatt (ID 2305)
hydroSites <- data.frame(x = c(723675, 718840, 727110, 737270),
                         y = c(252720, 248440, 247290, 251290))
widthObs <- c(35, 2.5, 20, 8) # estimated from aerial images
QObs <- c(19.96, 0.099, 3.195, 0.516) # mean values during 2012-2021

# locate hydrological stations to AG nodes. Use showPlot to visually inspect if attribution is correct
hydroAG <- numeric(length(hydroSites$x))
for (i in 1:length(hydroSites$x)){
  tmp <- rivnet::locate_site(hydroSites$x[i], hydroSites$y[i], Thur, showPlot = TRUE)
  title(i)
  hydroAG[i] <- tmp$AGnode
}
# site 3 is actually located upstream of the confluence: let's fix its attribution
upNodes <- which(Thur$AG$downNode==hydroAG[3])
hydroAG[3] <- upNodes[which.max(Thur$AG$A[upNodes])] # pick largest upstream reach

# assign hydraulic variables via hydro_river
hData <- c(widthObs, QObs)
hType <- c(rep("w", length(widthObs)), rep("Q", length(QObs)))
hNode <- c(rep(hydroAG,2))
Thur <- rivnet::hydro_river(data.frame(data = hData, type = hType, node = hNode),
                            Thur,
                            crossSection = "natural")

# WWTP data can be retrieved from https://map.geo.admin.ch/?lang=en&topic=ech&bgLayer=ch.swisstopo.pixelkarte-farbe&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bav.haltestellen-oev,ch.swisstopo.swisstlm3d-wanderwege,ch.astra.wanderland-sperrungen_umleitungen,ch.bafu.gewaesserschutz-klaeranlagen_ausbaugroesse&layers_opacity=1,1,1,0.8,0.8,0.75&layers_visibility=false,false,false,false,false,true&layers_timestamp=18641231,,,,,&catalogNodes=687,1743
ARAsites <- read.csv("ARAcoord.csv")
ARAnode <- numeric(length(ARAsites$n))
for (i in 1:length(ARAsites$n)){
  out <- rivnet::locate_site(ARAsites$X[i], ARAsites$Y[i], Thur, showPlot = TRUE)
  ARAnode[i] <- out$AGnode
  title(i)
  #readline(prompt="Press [enter] to continue:") # uncomment to inspect pictures one by one
}

# landcover data can be downloaded at https://www.bfs.admin.ch/bfs/de/home/dienstleistungen/geostat/geodaten-bundesstatistik/boden-nutzung-bedeckung-eignung/arealstatistik-schweiz.assetdetail.20104753.html
landCovCH <- read.csv("ag-b-00.03-37-area-csv.csv", sep = ";") # add this .csv file to the source folder
landCovCH <- terra::rast(data.frame(landCovCH$E, landCovCH$N, landCovCH$LU18_4),
                         type = "xyz", crs = "EPSG:2056") # convert into raster
# Legend: 1-Urban; 2-Agriculture; 3-Forest; 4-Improductive
landCovCH <- terra::project(landCovCH, crs("EPSG:21781")) # project into CH1903 coordinate system

# calculate covariate values via covariate_river
Thur <- rivnet::covariate_river(landCovCH, Thur)
names(Thur$SC$locCov) <- c("urban", "agriculture", "forest", "improductive")

# model parameters
k_diffuse <- c(0.3, 0.15, 0.1) # mgP m-2 day-1 (equiv to ~110, ~55, 36.5 mg/m2/yr for urban, agriculture, forested areas, respectively)
k_PopEq <- 150 # mg P/day/hab
vf <- 0.17 # m/day
phi_diffuse <- (k_diffuse[1]*Thur$SC$locCov$urban
                + k_diffuse[2]*Thur$SC$locCov$agriculture
                + k_diffuse[3]*Thur$SC$locCov$forest)*Thur$SC$A
phi_point <- numeric(Thur$AG$nNodes)
phi_point[ARAnode] <- ARAsites$PopEq*k_PopEq

# solve model in the downstream direction, one reach at a time
# (this is allowed because transport is only in the downstream direction)
sortedNodes <- sort(Thur$AG$AReach, index.return=T)
sortedNodes <- sortedNodes$ix

P_eq <- numeric(Thur$AG$nNodes)
for (i in 1:Thur$AG$nNodes){
  node <- sortedNodes[i]
  upNodes <- which(Thur$AG$downNode==node) # is node a headwater?
  if (length(upNodes)>0){
    upstreamP <- sum(86400*Thur$AG$discharge[upNodes]/Thur$AG$volume[upNodes]*P_eq[upNodes])
  } else {upstreamP <- 0}
  P_eq[node] <- (phi_diffuse[node] + phi_point[node] + upstreamP)/
    (86400*Thur$AG$discharge[node]/Thur$AG$volume[node] + vf/Thur$AG$depth[node])
}

# draw figure
pdf("Fig5.pdf", width = 21/2.54, height = 12/2.54)
par(mfrow = c(1,3))
kol <- colorRampPalette(c("chartreuse3", "darkgoldenrod1", "chocolate", "gray80", "gray90"))
plot(Thur, type = "elev2D", colPalette = kol(1000),
     drawRiver = TRUE, thrADraw = Thur$thrA,
     riverColor = "#0066FF",
     min_lwd = 0.2, max_lwd = 2)
points(hydroSites$x, hydroSites$y, pch = 19)
title ("Elevation [m a.s.l.]")

kols <- hcl.colors(1000, "viridis", rev = T)
log_phi_diffuse <- log10(phi_diffuse)
log_phi_diffuse[log_phi_diffuse==-Inf] <- NaN
plot(log_phi_diffuse, Thur, type = "SC",
     riverColor = "white", colPalette = kols,
     colLevels = c(4,7), min_lwd = 0.2, max_lwd = 2)
points(ARAsites$X, ARAsites$Y,
       bg = kols[round((log10(phi_point[ARAnode]) - 4)/3*1000)],
       pch = 21)
title("P loads [mg P /day]")

plot(P_eq/Thur$AG$volume, Thur, colLevels = c(0, 150),
     min_lwd = 0.2, max_lwd = 2)
title('P concentration [mg P/m3]')
dev.off()

