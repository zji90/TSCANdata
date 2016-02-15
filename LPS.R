library(monocle)
library(data.table)

###LPS 
rawdata <- fread("paper/rawdata/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt",header=T)
rn <- as.matrix(rawdata[,1,with=F])
rawdata <- as.matrix(rawdata[,-1,with=F])
row.names(rawdata) <- rn[,1]
rawdata <- rawdata[,!(log(rawdata["LYZ1",]+1) < 6 | log(rawdata["SERPINB6B",]+1) > 4)]
LPSdata1 <- as.matrix(rawdata[,grep("^LPS_1h_S[0-9]*$",colnames(rawdata))])
LPSdata2 <- as.matrix(rawdata[,grep("^LPS_2h_S[0-9]*$",colnames(rawdata))])
LPSdata4 <- as.matrix(rawdata[,grep("^LPS_4h_S[0-9]*$",colnames(rawdata))])
LPSdata6 <- as.matrix(rawdata[,grep("^LPS_6h_S[0-9]*$",colnames(rawdata))])
LPSdata <- cbind(LPSdata1,LPSdata2,LPSdata4,LPSdata6)
LPSdata <- log2(LPSdata + 1)
colnames(LPSdata) <- c(paste0("T1_",1:sum(grepl("^LPS_1h_S[0-9]*$",colnames(LPSdata)))),paste0("T2_",1:sum(grepl("^LPS_2h_S[0-9]*$",colnames(LPSdata)))),paste0("T3_",1:sum(grepl("^LPS_4h_S[0-9]*$",colnames(LPSdata)))),paste0("T4_",1:sum(grepl("^LPS_6h_S[0-9]*$",colnames(LPSdata)))))
data <- LPSdata[rowSums(LPSdata) > 0,]
clures <- hclust(dist(data))
cluster <- cutree(clures,0.05*nrow(data))
aggdata <- aggregate(data,list(cluster),mean)
aggdata <- aggdata[,-1]
#waterfall
LPSwaterfallorder <- pseudotimeprog.foo(aggdata,color="black") 

#kmeans
LPSkmeans <- exprkmeans(aggdata)
LPSkmeansorder <- TSCANorder(LPSkmeans)
