library(monocle)
library(data.table)

###HSMM rawdata
data(HSMM)
HSMMdata <- exprs(HSMM)
HSMMdata <- log2(HSMMdata+1)
colnames(HSMMdata) <- c(paste0("T1_",1:sum(grepl("T0",colnames(HSMMdata)))),paste0("T2_",1:sum(grepl("T24",colnames(HSMMdata)))),paste0("T3_",1:sum(grepl("T48",colnames(HSMMdata)))),paste0("T4_",1:sum(grepl("T72",colnames(HSMMdata)))))
save(HSMMdata,file="paper/HSMM/HSMMdata.rda")

###LPS rawdata
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
save(LPSdata,file="paper/LPS/LPSdata.rda")

###neural rawdata
neuraldata <- read.table("paper/rawdata/Single_TPM.txt",as.is=T,sep="\t",row.names = 1)
neuraldata <- as.matrix(log2(neuraldata+1))
save(neuraldata,file="paper/neural/neuraldata.rda")
