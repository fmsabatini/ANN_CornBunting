## ---- eval=F, echo=F-------------------------------------------------------------------------------------------
## #Start live server
## # Code for the instructor
## s <- livecode::serve_file(bitly=FALSE)
## system(paste0("ngrok http ", s$url), wait = FALSE, invisible = FALSE,
##        show.output.on.console = FALSE, minimized = FALSE)


## ---- message=F, warning=F-------------------------------------------------------------------------------------
# Install packages, if needed
#install.packages(c("dplyr", "readxl", "sf", "ggplot2", 
#                       "neuralnet", "corrplot", "cowplot"))

#Load Required libraries
# improved data manipulation
library(dplyr)
library(corrplot)
library(readxl)
# Spatial 
library(sf)
# visualization
library(ggplot2)
library(cowplot)
# neural networks
library(neuralnet)


## ----import data, cache=T, warning=F, message=F, results = 'hide'----------------------------------------------
## import data
Biotopes_sf <- st_read("data/Brandenburg_Biotops_2009.gpkg")
grid_sf <- st_read("data/zt_v100_UTM33N.gpkg")
mydata <- read_xls("data/zt_habitat_corn_bunting.xls",sheet = 1)


## ----check data1-----------------------------------------------------------------------------------------------
head(mydata)


## --------------------------------------------------------------------------------------------------------------
summary(mydata)
table(mydata$class)


## ----standardize predictors------------------------------------------------------------------------------------
## standardize variables to 0 - 1 range
# we drop the first (index) and the last two(HSI and response variable)columns
predictors <- mydata[,c(2:8)] 
maxs <- apply(predictors, 2, max) 
mins <- apply(predictors, 2, min)

scaled <- as.data.frame(scale(predictors, center = mins, scale = maxs - mins))
mydata.scaled <- data.frame(class=mydata$class, scaled)
summary(mydata.scaled)


## ----plot spatial data, cache=T--------------------------------------------------------------------------------
#plot data 
#Careful, the graph below renders slowly
(Landuse <- ggplot(data=Biotopes_sf) +  ## renders slowly
  geom_sf(aes(fill=Biotyp_8st)) + # use column biotyp_8st as a color code
  theme_classic() + 
  theme(legend.position = "none") #+ 
  #geom_sf(data=grid_sf, fill=NA)
)


## --------------------------------------------------------------------------------------------------------------
head(Biotopes_sf)


## ----split train vs test data----------------------------------------------------------------------------------
set.seed(899) # set a seed, so to make the resampling procedure reproducible
mydata.subset <- sample_frac(mydata.scaled, size=0.2) 
#function sample_frac comes from dplyr package! 

#split in train and test
n <- nrow(mydata.subset)
share.train <- 0.8
train.id <- sample(1:n, n*share.train, replace=F) 
train = mydata.subset[train.id,]
test = mydata.subset[-train.id,]


## --------------------------------------------------------------------------------------------------------------
myformula <- as.formula("class ~ NS + TM + WN + SN + AK + DL + WW")


## ----train network (0-1), cache=T------------------------------------------------------------------------------
set.seed(900)
nn_5 <- neuralnet(myformula, train, hidden=c(5), lifesign='minimal', threshold=0.01,
                  linear.output = F)


## ----train nn on 0.8-0.2, cache=T------------------------------------------------------------------------------
train_recoded <- train
train_recoded$class <- train_recoded$class*0.6+0.2
table(train_recoded$class)
nn_5 <- neuralnet(myformula, train_recoded, hidden=c(5), lifesign="full", threshold=0.01,
                  linear.output = F)


## ----plot nn output, fig.align="center", fig.height=10, fig.width=10, message=F, warning=F---------------------
plot(nn_5, rep = "best")


## ----check weights---------------------------------------------------------------------------------------------
nn_5$weights


## ----Partial Dependency Plots, fig.align="center", fig.height=3, fig.width=4, message=F, warning=F-------------
# create a new dataset
newdata <- as.data.frame(t(apply(scaled, MARGIN=2, "mean")))
newdata <- newdata[rep(1,100),] #repeat the means 100 times
newdata$AK <- seq(0,1, length.out = 100)
# feed the new dataset to the NN, and predict the outputs
predict_newdata = neuralnet::compute(nn_5, newdata)
newdata$predict = predict_newdata$net.result
# plot
ggplot(data=newdata) + 
  geom_line(aes(x=AK, y=predict)) + 
  theme_bw() 


## ----calculate MSE, fig.align="center", fig.height=4, fig.width=5, message=F, warning=F------------------------
#transform test data to 0.2 - 0.8
test_recoded <- test
test_recoded$class <- test_recoded$class*0.6+0.2

## Feed the test data into our ANN, and 'predict' the output
predict_testNN = neuralnet::compute(nn_5, test_recoded)
predict_testNN = predict_testNN$net.result 

# Visual comparison of predicted vs expected values
plot(test$class, predict_testNN, col='blue', pch=16, ylab = 
       "predicted class NN", xlab = "real class")
#same in tabular data
predict_01 <- ifelse(predict_testNN<0.5, 0, 1)
table(predict_01, test$class)

# Calculate Root Mean Square Error (RMSE)
SSE = sum((test_recoded$class - predict_testNN)^2)
RMSE = (SSE / nrow(test)) ^ 0.5


## ----function for CV-------------------------------------------------------------------------------------------
mydata.subset$class <- mydata.subset$class*0.6+0.2 #recode to improve performance

get.RMSE <- function(i, data, form, hidden, share.train=0.8, label="full"){
  train.id <- sample(1:nrow(data), n*share.train, replace=F) 
  train = data[train.id,]
  test = data[-train.id,]
  nn <- neuralnet(form=form, train, hidden=c(hidden), 
                  lifesign = "full", threshold=0.01, linear.output = F)
  predict_testNN <- neuralnet::compute(nn, test)
  predict_testNN = predict_testNN$net.result 
  RMSE.NN = (sum((test$class - predict_testNN)^2) / nrow(test)) ^ 0.5
  #print(i)
  return(data.frame(config=paste(label, paste(hidden, collapse="_")), RMSE=RMSE.NN))
}


## ----test different configurations, eval=F---------------------------------------------------------------------
## #runs in about ~10-15 minutes
## set.seed(200)
## RMSE.2 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(2))
## RMSE.3 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(3))
## RMSE.4 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(4))
## RMSE.5 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(5))
## RMSE.6 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(6))
## 
## RMSE.2.2 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(2,2))
## RMSE.3.3 <- lapply(1:5, get.RMSE, mydata.subset, myformula, c(3,3))
## ## compile and visualize output
## RMSE.df <- do.call(rbind, c(RMSE.2, RMSE.3, RMSE.4, RMSE.5, RMSE.6, RMSE.2.2, RMSE.3.3))


## ---- echo=F---------------------------------------------------------------------------------------------------
load("RMSE.Rdata")


## ----Correlation matrix, fig.align="center", fig.height=4, fig.width=5, message=F, warning=F-------------------
res <- cor(scaled, use = "pairwise.complete.obs")
corrplot(res, type = "upper", 
         tl.col = "black", tl.srt = 45, number.cex=0.6, addCoef.col = "black", diag=F)



## ----reduced NN, eval=F----------------------------------------------------------------------------------------
## 
## RMSE.3.no_AK <- lapply(1:5, get.RMSE, mydata.subset,
##                        as.formula("class ~ NS + TM + WN + SN + DL + WW"), c(3),
##                        label="no_AK")
## RMSE.3.no_DL <- lapply(1:5, get.RMSE, mydata.subset,
##                        as.formula("class ~ NS + TM + WN + SN + AK + WW"), c(3),
##                        label="no_DL")
## RMSE.df2 <- rbind(RMSE.df, do.call(rbind, c(RMSE.3.no_AK, RMSE.3.no_DL)))


## ----plot CV output, fig.align="center", fig.height=4, fig.width=5, message=F, warning=F-----------------------
ggplot(data=RMSE.df2) + 
  geom_boxplot(aes(x=config, y=RMSE))


## ----Predict over Study area-----------------------------------------------------------------------------------
predict_alldata = neuralnet::compute(nn_5, mydata.scaled)
predict_01 <- ifelse(predict_alldata$net.result<0.5, 0, 1)
predict_01 <- data.frame(PRED01=predict_01, 
                         ZT_V100_=mydata$ZT_V100_ID)
# join output to fishnet grid
grid_sf2 <- left_join(grid_sf, predict_01, by="ZT_V100_")


## ----Predict over study area reduced NN------------------------------------------------------------------------
myformula2 <- class ~ NS + TM + WN + SN + DL + WW
nn_3_noAK <- neuralnet(form=myformula2, train_recoded, hidden=c(3), 
                       lifesign = "full", threshold=0.01, linear.output = F)
predict_alldata = neuralnet::compute(nn_3_noAK, mydata.scaled)
predict_01_noAK <- ifelse(predict_alldata$net.result<0.5, 0, 1)
predict_01_noAK <- data.frame(PRED01=predict_01_noAK, 
                         ZT_V100_=mydata$ZT_V100_ID)

grid_sf3 <- left_join(grid_sf, predict_01_noAK, by="ZT_V100_")


## ----Plotting--------------------------------------------------------------------------------------------------
HSI <- mydata[,c(1,9)]
names(HSI)[1] <- "ZT_V100_"
grid_sf_hsi <- left_join(grid_sf, HSI,  by="ZT_V100_")


## --------------------------------------------------------------------------------------------------------------
pred.full <- ggplot(data=grid_sf2) +
  geom_sf(aes(fill=as.factor(PRED01)), alpha=0.3, col=NA) + 
  labs(fill='Suitable') +
  theme_bw() + 
  theme(axis.text = element_blank())

pred.noAK <- pred.full %+% grid_sf3 #update graph with new data

HSI <- ggplot(data=grid_sf_hsi) + 
  geom_sf(aes(fill=HSI), col=NA) +
  theme_bw() + 
  theme(axis.text = element_blank())


## ---- fig.align="center", fig.height=8, fig.width=8, message=F, warning=F, cache=T-----------------------------
cowplot::plot_grid(pred.full, pred.noAK,
                   HSI, Landuse + theme(axis.text = element_blank()), 
                   nrow=2, rel_heights = c(1,1,1,0.6))


## ----Overlap prediction with landuse map, fig.align="center", fig.height=4, fig.width=8, message=F, warning=F, cache=T----
left <- ggplot(data=Biotopes_sf) + 
  geom_sf(fill=NA) + 
  geom_sf(data=grid_sf2, aes(fill=as.factor(PRED01)), alpha=0.5, col=NA) + 
  labs(fill='Full ANN') +
  theme_bw() + 
  theme(axis.text = element_blank())
right <- ggplot(data=Biotopes_sf) + 
  geom_sf(fill=NA) + 
  geom_sf(data=grid_sf3, aes(fill=as.factor(PRED01)), alpha=0.5, col=NA) + 
  labs(fill='Red. ANN') +
  theme_bw() + 
  theme(axis.text = element_blank())

cowplot::plot_grid(left, right, nrow=1)




## ----save output, warning=F, message=F, eval=T, results="hide"-------------------------------------------------
dir.create("_output")
st_write(grid_sf2, dsn="_output/CornBunting_NN3_pred.shp", append=T)
st_write(grid_sf3, dsn="_output/CornBunting_NN3_noAK_pred.shp", append=T)

write.csv(predict_01, file = "_output/CornBunting_NN3_pred01.csv")
write.csv(predict_01_noAK, file = "_output/CornBunting_NN3_noAK_pred01.csv")


## --------------------------------------------------------------------------------------------------------------
sessionInfo()

