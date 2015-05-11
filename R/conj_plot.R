rm(list=ls())
gc()

library(Rcpp)
library(RcppEigen)
library(bayesGDS)
library(reshape2)
library(ggplot2)
library(gdata)

facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)

  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))

  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }

  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}


mod.name <- "conjugate"
theme_set(theme_bw())

file.path <- paste0("~/Documents/gds-mksc/",mod.name,"/")

mode.file <- paste0(file.path,"results/",mod.name,"_mode.Rdata")
draws.file <- paste0(file.path,"results/",mod.name,"_all.Rdata")
boxplot.file <- "~/Documents/gds-mksc/paper/rev1/conjugate_box.pdf"


load(mode.file)
load(draws.file)

tmp <- dcast(draws,iterations+method~parameters)
tmp <- transform(tmp,tau2=exp(2*log.tau),sig2=exp(2*log.sig))
draws <- melt(tmp, id.vars=c("iterations","method"),
              variable.name="parameters")


vars <- c("mu","sig2","tau2","log.post")

d <- subset(draws, parameters %in% vars & iterations>0)
d$parameters <- reorder(d$parameters,new.order=c("mu","tau2","sig2","log.post"))

d$method <- mapvalues(d$method, from=c("mcmc","gds"), to=c("Gibbs","Ours"))


P <- ggplot(data=d,aes(x=value, linetype=method)) + geom_density()
P <- P + facet_wrap(~parameters, scales="free",nrow=1,ncol=4, as.table=FALSE)

P <- P + theme(legend.position="top",
               axis.text.y=element_blank(),
               axis.title.y=element_blank(),
               
               axis.text=element_text(size=rel(.6)))
P <- facet_wrap_labeller(P,labels=c(expression(mu),
                               expression(tau^2),
                               expression(sigma^2),
                               expression("log posterior"))
                         )
##print(P)                                

B <- ggplot(data=d,aes(x=method,y=value)) + geom_boxplot(outlier.size=.6,size=.2)
B <- B + theme(axis.text=element_text(size=8),
##               axis.title=element_text(size=10),
               axis.title=element_blank(),
               strip.background=element_rect(fill="white"))
B <- B + facet_wrap(~parameters, scales="free",nrow=1,ncol=4)
B <- facet_wrap_labeller(B,labels=c(expression(mu),
                               expression(tau^2),
                               expression(sigma^2),
                               expression("log posterior"))
                         )


pdf(file=boxplot.file, width=6,height=2)
print(B)
dev.off()

## trace <- ggplot(d,aes(x=iterations,y=value)) + geom_line()
## trace <- trace + facet_wrap(~parameters, scales="free")
##print(trace) 


