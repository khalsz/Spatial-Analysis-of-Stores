# Spatial Econometrics Project #
# M1 DASEE students: Khalid, Siarhei, Shahkar, Sakina

library(spdep)
library(rgdal)
library(tidyverse)
library(sf)

library(classInt)
library(maps)         
library(maptools)
library(spgwr)
library(grid)
library(gridExtra)

library(ggplot2)

library(sp)
library(rgeos)
library(tmap)
library(tmaptools)
library(spgwr)
library(grid)
library(gridExtra)

     ## Data management
library(sp)           ## Data management
library(spdep)        ## Spatial autocorrelation
library(gstat)        ## Geostatistics
library(splancs)      ## Kernel Density
library(spatstat)     ## Geostatistics
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(spgwr)
library(spatialreg)

setwd ("C:/Users/Siarhei/Downloads")
Rdvf <- read.csv("C:/Users/Siarhei/Downloads/siren.dvf.final.csv")
summary(Rdvf)

#extracting coordinate
coordinates(Rdvf) <- ~latitude+longitude
cor <- coordinates(Rdvf)

#converting into point shapefile
my.sf.poin <- st_as_sf(x = Rdvf, coords = c("cor"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# matrix with the indices of points belonging to the set of the k nearest neighbors
k1 <- knn2nb(knearneigh(cor, k = 1))

# calculating upper bound distance
Eucl_k1dists<- max(unlist(nbdists(k1,cor)))

#Euclidean Distance neigbours 
nb.dist.band <- dnearneigh(cor, 0, Eucl_k1dists)
nb.dist.band

# Euclidean distance between linked regions
distances <- nbdists(nb.dist.band,cor)

# inverse of the Euclidean distance
invd1 <- lapply(distances, function(x) (1/x))

# Spatial weight matrix of Inverse distance neigbours
invd.weights <- nb2listw(nb.dist.band, glist = invd1, style = "W", zero.policy = TRUE)
invd.weights

# Spatial weight matrix of Euclidean distance neigbours
rswm_q <- nb2listw(nb.dist.band, style = "W", zero.policy = TRUE)
rswm_q

#Plot showing links between Neighbors (euclidean distance)
plot(nb.dist.band, cor, lwd=.2, col="blue", cex = .5)
#Plot showing links between Neighbors (inverse euclidean distance)
plot(invd.weights, cor, lwd=.2, col="red", cex = .5)
 
# card function tallies the number of neighbors.
dist.band.card <- card(nb.dist.band)

ggplot() +
  geom_bar(aes(x=dist.band.card)) +
  xlab("Number of Neighbors")

#to change scientific format
options(scipen = 9)
#Moran for Euclidean distance weight
moran.test(my.sf.poin$value, listw = rswm_q, zero.policy = TRUE, na.action = na.omit)

#Monte-Carlo Moran for Euclidean distance weight
set.seed(1234)
Ebperm = moran.mc(my.sf.poin$value, listw = rswm_q, nsim = 999, zero.policy = TRUE, na.action = na.omit)
Ebperm

#Hitsogram chart (Monte-Carlo Moran for Euclidean distance weight)
hist(Ebperm$res, freq = TRUE, breaks = 20, xlab = "Simulated Moran's I")
abline(v=0, col="red")
plot(Ebperm, main="", las=1)

#Moran for inverse distance weight
moran.test(my.sf.poin$value, listw = invd.weights, zero.policy = TRUE, na.action = na.omit)

#Monte-Carlo Moran for inverse distance weight
set.seed(1234)
bperm = moran.mc(my.sf.poin$value, listw = invd.weights, nsim = 999, zero.policy = TRUE, na.action = na.omit)
bperm

#Hitsogram chart (Monte-Carlo Moran for inverse distance weight)
hist(bperm$res, freq = TRUE, breaks = 20, xlab = "Simulated Moran's I")
abline(v=0, col="red")

plot(bperm, main="", las=1)

#Moran Scatter Plot for Euclidean Distance 
moran <- moran.plot(my.sf.poin$value, listw = rswm_q, main = "Moran Scatterplot for Euclidean Distance")

#Moran Scatter Plot for Inverse Distance 
moran <- moran.plot(my.sf.poin$value, listw = invd.weights, main = "Moran Scatterplot for Inverse Euclidean Distance")

################ Geary C ################ 

#Geary C for Euclidean distance weight
geary.test(my.sf.poin$value, listw = rswm_q)

#Monte-Carlo Geary C for Euclidean distance weight
set.seed(1234)
Ebper = geary.mc(my.sf.poin$value, listw = rswm_q, nsim = 999)
Ebper

#Hitsogram chart (Monte-Carlo Geary C for Euclidean distance weight)
hist(Ebper$res, freq = TRUE, breaks = 20, xlab = "Simulated Geary's C")
abline(v=0, col="red")

plot(Ebper, main="", las=1)

#Geary C for inverse distance weight
geary.test(my.sf.poin$value, listw = invd.weights)

#Monte-Carlo Geary C for inverse distance weight
set.seed(1234)
bper = geary.mc(my.sf.poin$value, listw = invd.weights, nsim = 999)
bper

#Hitsogram chart (Monte-Carlo Geary C for inverse distance weight)
hist(bper$res, freq = TRUE, breaks = 20, xlab = "Simulated Geary's C")
abline(v=0, col="red")
plot(bper, main="", las=1)

############ Local Moran ##########

#Local Moran for Euclidean Distance Matrix
Elocal <- localmoran(x = my.sf.poin$value, listw = rswm_q)
summary(Elocal)
Elocl <- as.data.frame(Elocal)

# Map the local Moran's I values (Euclidean Distance)
Efj5 <- classIntervals(round(Elocl$Ii, digits = 4), n = 5, style = "fisher")
Epal <- grey.colors(4, 0.95, 0.55, 2.2)
Efj5Colours <- findColours(Efj5, Epal)
plot(Elocal, col = Efj5Colours, pch = 19, border = NA)
legend("bottomright", fill = attr(Efj5Colours, "palette"), legend = names(attr(Efj5Colours, 
                                                                              "table")), bty = "n")

#Local Moran for Inverse Distance Matrix
local <- localmoran(x = my.sf.poin$value, listw = invd.weights)
summary(local)
locl <- as.data.frame(local)

# Map the local Moran's I values (Inverse Distance)
fj5 <- classIntervals(round(locl$Ii, digits = 4), n = 5, style = "fisher")
pal <- grey.colors(4, 0.95, 0.55, 2.2)
fj5Colours <- findColours(fj5, pal)
plot(local, col = fj5Colours, pch = 19, border = NA)
legend("bottomright", fill = attr(fj5Colours, "palette"), legend = names(attr(fj5Colours, 
                                                                              "table")), bty = "n")

#getis ord for inverse distance
local_g <- localG(my.sf.poin$value, invd.weights)
local_g <- cbind(my.sf.poin, as.matrix(local_g))
head(local_g)
names(local_g)[6] <- "gstat"

getis = as.vector(localG(my.sf.poin$value, invd.weights))

############ Spatial regression ############ 

#converting categorical variable into factor
Rdvf$group = as.factor(Rdvf$group)
summary(Rdvf$group)

#OLS regression
reqeq1 = value ~ group
options(scipen = 9)

reg1 = lm(reqeq1, data = Rdvf)
summary(reg1)

#moran test 
lm.morantest(reg1, rswm_q)

#Lagrange Multiplier Test
lm.LMtests(reg1, rswm_q, test = "all")

#lag model has the lowest p-value, so it is the best SLX. Also, it has the highest spillover effect
#i.e RLMlag = 8.6495, df = 1, p-value = 0.003272

#SLX Spatially lagged X, Equation: Y (shop price) = X(shop activity)B (beta)+ WX (weight matrix of x) T (tetha) +e
reg2 = lmSLX(reqeq1, data = Rdvf, rswm_q )
summary(reg2)

# Result of direct and indirect and its significance
summary (impacts(reg2, listw = rswm_q, R = 500), zstats=TRUE)
