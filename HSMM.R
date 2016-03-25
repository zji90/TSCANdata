library(monocle)
library(data.table)

###HSMM without prior information
data(HSMM)
HSMMdata <- exprs(HSMM)
HSMMdata <- log2(HSMMdata+1)
colnames(HSMMdata) <- c(paste0("T1_",1:sum(grepl("T0",colnames(HSMMdata)))),paste0("T2_",1:sum(grepl("T24",colnames(HSMMdata)))),paste0("T3_",1:sum(grepl("T48",colnames(HSMMdata)))),paste0("T4_",1:sum(grepl("T72",colnames(HSMMdata)))))
data <- HSMMdata
data <- data[rowSums(data) > 0,]
clures <- hclust(dist(data))
cluster <- cutree(clures,0.05*nrow(data))
aggdata <- aggregate(data,list(cluster),mean)
aggdata <- aggdata[,-1]
#waterfall
HSMMwaterfallorder <- pseudotimeprog.foo(aggdata,color="black") 

#kmeans
HSMMkmeans <- exprkmeans(aggdata)
HSMMkmeansorder <- TSCANorder(HSMMkmeans)

#Scuba
fitpc <- principal.curve(t(aggdata), smoother = "lowess", maxit = 200)
HSMMscubaorder <- colnames(aggdata)[order(fitpc$lambda)]

