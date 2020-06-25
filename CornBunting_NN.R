library(dplyr)
library(raster)
library(rgdal)
library(sf)
library(ggplot2)
library(neuralnet)


## import data
Biotopes <- readOGR("Brandenburg_Biotops_2009.gpkg")#, encoding = "UTF-8", use_iconv=T) # encoding not behaving as expected
mydata <- readxl::read_xls("zt_habitat_corn_bunting.xls",sheet = 1)
grid <- readOGR("zt_v100_UTM33N.gpkg")

### Explore input data 
## Quick visualization of the data
grid_sf <- st_as_sf(grid)
Biotopes_sf <- st_as_sf(Biotopes)
ggplot(data=Biotopes_sf) +  ## renders slowly
  geom_sf(aes(fill=Biotyp_8st)) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  geom_sf(data=grid_sf, fill=NA)

## check data ranges and type
summary(mydata)

## standardize variables to -1 1 range
maxs <- apply(select(mydata, -ZT_V100_ID, -class), 2, max) 
mins <- apply(select(mydata, -ZT_V100_ID, -class), 2, min)

scaled <- as.data.frame(scale(select(mydata, -ZT_V100_ID, -class), 
                              center = mins, scale = maxs - mins))
mydata.scaled <- data.frame(class=mydata$class, scaled)
mydata.scaled.08 <- mydata.scaled %>% 
  mutate(class=ifelse(class==1, 0.8, 0.2))

# Split data in training (80%) and test (20%) subsets
set.seed(899)
n <- nrow(mydata.scaled)
train.id <- sample(1:n, n*.1, replace=F)  ### set to 0.1 to reduce computing time
train = mydata.scaled[train.id,]
train.08 = mydata.scaled.08[train.id,]
test = mydata.scaled[-train.id,]
test.08 = mydata.scaled.08[-train.id,]

myformula.char <- "class ~ NS + TM + WN + SN + AK + DL + WW"
myformula <- as.formula(myformula.char)
nn_5 <- neuralnet(myformula, train, hidden=c(5), lifesign="full", threshold=0.01, linear.output = F)
nn_5.08 <- neuralnet(myformula, train.08, hidden=c(5), lifesign="full", threshold=0.01, linear.output = F)


nn_3 <- neuralnet(myformula, train, hidden=c(3), lifesign="full", threshold=0.1, linear.output = F)
nn_5_3 <- neuralnet(myformula, train, hidden=c(5,3), lifesign="full", threshold=0.1, linear.output = F)
nn_5_5 <- neuralnet(myformula, train, hidden=c(5,3), lifesign="full", threshold=0.1, linear.output = F)

plot(nn_5)
plot(nn_5_3)


train.list <-  lapply(1:10, function(x){sample(1:n, n*.8, replace=F)})

get.RMSE <- function(data, train.index, form, hidden){
  train.data <- data[train.index,]
  test.data <- data[-train.index,]
  nn <- neuralnet(form, train, hidden=c(hidden), lifesign="full", threshold=0.01, linear.output = F)
  predict_testNN <- neuralnet::compute(nn, test)
  predict_testNN = (predict_testNN$net.result * (max(data$class) - min(data$class))) + min(data$class)
  RMSE.NN = (sum((test$class - predict_testNN)^2) / nrow(test.data)) ^ 0.5
  return(RMSE.NN)
}

rmse_h3 <- get.RMSE(mydata.scaled, train.list[[1]], form=myformula, hidden=c(3))
rmse_h3_list <- lapply(train.list, get.RMSE, data=mydata.scaled, form=myformula, hidden=c(3))
rmse_h5_list <- lapply(train.list, get.RMSE, data=mydata.scaled, form=myformula, hidden=c(5))






## Prediction using neural network

predict_testNN = neuralnet::compute(nn_5, test)
predict_testNN = (predict_testNN$net.result * (max(mydata.scaled$class) - min(mydata.scaled$class))) + 
  min(mydata.scaled$class)

plot(test$class, predict_testNN, col='blue', pch=16, ylab = 
       "predicted class NN", xlab = "real class")
#abline(0,1)

predict_01 <- ifelse(predict_testNN<0.5, 0, 1)
table(predict_01, test$class)

# Calculate Root Mean Square Error (RMSE)
RMSE.NN = (sum((test$class - predict_testNN)^2) / nrow(test)) ^ 0.5


  

