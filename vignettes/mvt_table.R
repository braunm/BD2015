rm(list=ls())
gc()

library(plyr)
library(reshape2)
library(ggplot2)
library(xtable)

load.file <- "~/Documents/Direct Sampling/results/rev1/MVT/LL.Rdata"
save.file <- "~/Documents/Direct Sampling/paper/rev1/LL.tex"

load(load.file)

cols <- c("iter","n.x","n.obs","scale","M","LML.ds","log.hme","LML.true","LML.Lenk","log.hme","acc.rate","time")

w <- which(!is.na(res[,1]))

R <- data.frame(res[w,cols])
R <- transform(R, mape.ds=100*abs(LML.ds-LML.true)/abs(LML.true),
               mape.lenk=100*abs(LML.Lenk-LML.true)/abs(LML.true),
               acc.pct=100*acc.rate)

fr <- melt(R,id.vars=c("iter","n.x","n.obs","scale","M"))

tab.mean <- acast(fr, n.x+n.obs+M+scale~variable, mean, na.rm=TRUE)
tab.sd <- acast(fr, n.x+n.obs+M+scale~variable, sd,na.rm=TRUE)

ids <- colsplit(rownames(tab.mean),"_",paste("V",1:4,sep=""))

tab.mean <- cbind(ids,tab.mean)
tab.sd <- cbind(ids,tab.sd)

tab1 <- cbind(ids,
              tab.mean[,"LML.true"],tab.sd[,"LML.true"],
              tab.mean[,"LML.ds"],tab.sd[,"LML.ds"],
              tab.mean[,"LML.Lenk"],tab.sd[,"LML.Lenk"],
              tab.mean[,"log.hme"],tab.sd[,"log.hme"]
              )

tab2 <- cbind(ids,
              tab.mean[,"mape.ds"],tab.sd[,"mape.ds"],
              tab.mean[,"mape.lenk"],tab.sd[,"mape.lenk"],
              tab.mean[,"acc.pct"],tab.sd[,"acc.pct"],
              tab.mean[,"time"], tab.sd[,"time"]
              )
 

colnames(tab1) <- c("k","n","M","scale",
                    "mean","sd",
                    "mean","sd",
                    "mean","sd",
                    "mean","sd"
                    )


colnames(tab2) <- c("k","n","M","scale",
                    "mean","sd",
                    "mean","sd",
                    "mean","sd",
                    "mean","sd"
                    )

x1 <- xtable(tab1,
             align="l|llll|rrrrrrrr|",
             digits=c(0,1,1,1,1,0,1,0,1,0,1,0,1),
             display=c("s",
               "d","d","d","f",
               rep("f",8)
             )
             )


x2 <- xtable(tab2,
             align="l|llll|rrrrrrrr|",
             digits=c(0,1,1,1,1,2,2,2,2,2,2,2,2),
             display=c("s",
               "d","d","d","f",
               rep("f",8)
               )
             )
             
print(x1, file=save.file,
      include.rownames=FALSE,
      hline.after=c(-1,0,seq(4,nrow(tab1),by=4),nrow(tab1)),
      floating=TRUE,
      floating.environment="table"
      )


print(x2, file=save.file,
      include.rownames=FALSE,
      hline.after=c(-1,0,seq(4,nrow(tab2),by=4),nrow(tab2)),
      floating=TRUE,
      floating.environment="table",
      append=TRUE
      )

             

             
      
    

