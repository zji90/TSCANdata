library(monocle)
library(data.table)

###qNSC 
neuraldata <- read.table("paper/rawdata/Single_TPM.txt",as.is=T,sep="\t",row.names = 1)
neuraldata <- as.matrix(log2(neuraldata+1))
data <- neuraldata[rowSums(neuraldata) > 0,]
clures <- hclust(dist(data))
cluster <- cutree(clures,0.05*nrow(data))
aggdata <- aggregate(data,list(cluster),mean)
aggdata <- aggdata[,-1]
#waterfall
neuralwaterfallorder <- pseudotimeprog.foo(aggdata,color="black") 

#kmeans
neuralkmeans <- exprkmeans(aggdata)
neuralkmeansorder <- TSCANorder(neuralkmeans)

#Scuba
fitpc <- principal.curve(t(aggdata), smoother = "lowess", maxit = 200)
neuralscubaorder <- colnames(aggdata)[order(fitpc$lambda)]
