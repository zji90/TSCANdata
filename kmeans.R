exprkmeans <- function(data, clusternum = NULL, reduce = T) { 
      if (reduce) {
            sdev <- prcomp(t(data),scale=T)$sdev[1:20]
            x <- 1:20
            optpoint <- which.min(sapply(2:19, function(i) {
                  x2 <- pmax(0,x-i)
                  sum(lm(sdev~x+x2)$residuals^2)
            }))
            pcadim = optpoint + 1
            
            tmpdata <- t(apply(data,1,scale))
            colnames(tmpdata) <- colnames(data)
            tmppc <- prcomp(t(tmpdata),scale=T)
            pcareduceres <- t(tmpdata) %*% tmppc$rotation[,1:pcadim]      
      } else {
            pcareduceres <- t(data)
      }
      
      if (is.null(clusternum)) {
            ssper <- sapply(1:20, function(i) {
                  mean(sapply(1:5,function(time) {
                        res <- kmeans(pcareduceres,i)
                        res$betweenss/res$totss      
                  }))
            }) 
            
            x <- 1:20
            clusternum <- which.min(sapply(1:20, function(i) {
                  x2 <- pmax(0,x-i)
                  sum(lm(ssper~x+x2)$residuals^2)
            }))
      }
      
      res <- kmeans(pcareduceres,clusternum)
      clusterid <- res$cluster            
      clucenter <- res$centers      
      dp <- as.matrix(dist(clucenter))                                            
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      list(pcareduceres=pcareduceres,MSTtree=dp_mst,clusterid=clusterid,clucenter=clucenter)      
}
