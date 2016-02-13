library(monocle)
library(data.table)

###HSMM with prior information
data(HSMM)
HSMMdata <- exprs(HSMM)
HSMMdata <- log2(HSMMdata+1)
colnames(HSMMdata) <- c(paste0("T1_",1:sum(grepl("T0",colnames(HSMMdata)))),paste0("T2_",1:sum(grepl("T24",colnames(HSMMdata)))),paste0("T3_",1:sum(grepl("T48",colnames(HSMMdata)))),paste0("T4_",1:sum(grepl("T72",colnames(HSMMdata)))))
data <- HSMMdata
aggdata <- data[fData(HSMM)$use_for_ordering,]
HSMMwaterfallorder <- pseudotimeprog.foo(aggdata,color="black") 
