rm(list=ls())
gc()

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)

load.file <- "~/Documents/Direct Sampling/results/rev1/MVT/LL.Rdata"
save.file <- "~/Documents/Direct Sampling/paper/rev1/LL-short.tex"

load(load.file)

cols <- c("iter","n.x","n.obs","scale","M","LML.ds","log.hme","LML.true","LML.Lenk","log.hme","acc.rate","time")
w <- which(!is.na(res[,1]))

R <- data.frame(res[w,cols])
R <- transform(R, mape.ds=100*abs(LML.ds-LML.true)/abs(LML.true),
               mape.lenk=100*abs(LML.Lenk-LML.true)/abs(LML.true),
               acc.pct=100*acc.rate,
               secs=time*60)

fr <- melt(R,id.vars=c("iter","n.x","n.obs","scale","M"))
wh <- fr$scale %in% c(0.7,0.8) & fr$M==10000
fr <- fr[wh,]

tab.mean <- acast(fr, scale+n.x+n.obs~variable, mean, na.rm=TRUE)

ids <- colsplit(rownames(tab.mean),"_",paste("V",1:3,sep=""))

tab.mean <- cbind(ids,tab.mean)


tab1 <- cbind(ids,
              tab.mean[,"LML.true"],
              tab.mean[,"LML.ds"],
              tab.mean[,"log.hme"],  
              tab.mean[,"acc.pct"],
              tab.mean[,"secs"]
              
              )



colnames(tab1) <- c("scale","k","n",
                    "MVT",
                    "GDS",
                    "HME",
                    "acc %",
                    "secs"
                    )


x1 <- xtable(tab1,
             align="l|lll|rrrrr|",
             digits=c(0,1,1,0,0,0,0,1,1),
             display=c("s",
               "f","d","f",
               rep("f",5)
               )
             )


             
print(x1, file=save.file,
      include.rownames=FALSE,
      floating=TRUE,
      floating.environment="table"
      )



             
      
    

