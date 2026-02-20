## load libraries
library(phytools)
## read tree from file
anolis.tree<-read.tree("Anolis.tre")
print(anolis.tree,printlen=2)

anolis.data<-read.csv("anole.data.csv",row.names=1,stringsAsFactors=TRUE)
head(anolis.data)

## extract total body length and log-transform
lnSVL<-setNames(log(anolis.data$SVL),rownames(anolis.data))
head(lnSVL)


lnSVL[lnSVL>1.4] <- "big"
lnSVL[lnSVL<1.4] <- "small"
lnSVL <- factor(lnSVL, levels = c("small", "big"))
## set colors for plotting
cols<-setNames(c("red","lightblue"),levels(lnSVL))
## plot the tree & data
plotTree.datamatrix(anolis.tree,as.data.frame(lnSVL),colors=list(cols),header=FALSE,fsize=0.45)
## add legend
legend("topright",legend=levels(lnSVL),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

## fit ER model
fitER<-fitMk(anolis.tree,lnSVL,model="ER")
## fit ARD model
fitARD<-fitMk(anolis.tree,lnSVL,model="ARD")
## fit bite->suction model
fit01<-fitMk(anolis.tree,lnSVL,model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
## fit suction->bite model
fit10<-fitMk(anolis.tree,lnSVL,model=matrix(c(0,0,1,0),2,2,byrow=TRUE))
## extract AIC values for each model
aic<-c(AIC(fitER),AIC(fitARD),AIC(fit01),AIC(fit10))
## print summary table
data.frame(model=c("ER","ARD","small->big","big->small"),logL=c(logLik(fitER),logLik(fitARD),logLik(fit01),logLik(fit10)),AIC=aic,delta.AIC=aic-min(aic))

library(corHMM)
## create new data frame for corHMM
anolis.data<-data.frame(species=names(lnSVL),lnSVL=as.numeric(lnSVL))
head(anolis.data,n=10)
  ## estimate marginal ancestral states under a customed model
fit.marginal<-corHMM(anolis.tree,anolis.data,node.states="marginal",rate.cat=1,model = "ARD")
head(fit.marginal$states)

## plot the tree & data
plotTree.datamatrix(anolis.tree,as.data.frame(lnSVL),colors=list(cols),header=FALSE,fsize=0.45)
## add legend
legend("topright",legend=levels(lnSVL),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)
## add node labels showing marginal ancestral states
nodelabels(pie=fit.marginal$states,piecol=cols,cex=0.5)


## generate one stochastic character history
mtree<-make.simmap(anolis.tree,lnSVL,model="ARD")

## plot single stochastic map
plot(mtree,cols,fsize=0.4,ftype="i",lwd=2,offset=0.4,ylim=c(-1,Ntip(anolis.tree)))
## add legend
legend("bottomleft",legend=levels(lnSVL),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

## generate 1,000 stochastic character maps in which
## the transition rate is sampled from its posterior
## distribution
mtrees<-make.simmap(anolis.tree,lnSVL,model="ARD",nsim=100,Q="mcmc",vQ=0.01,prior=list(use.empirical=TRUE),samplefreq=10)
mtrees

## set plot margins
par(mar=c(5.1,4.1,2.1,2.1))
## create a plot of the posterior density from stochastic
## mapping
plot(d<-density(sapply(mtrees,function(x) x$Q[1,2]),bw=0.005),bty="n",main="",xlab="q",xlim=c(0,0.5),ylab="Posterior density from MCMC",las=1,cex.axis=0.8)
polygon(d,col=make.transparent("blue",0.25))
## add line indicating ML solution for the same parameter
abline(v=fit.marginal$solution[1,2])
text(x=fit.marginal$solution[1,2],y=max(d$y),"MLE(q)",pos=4)

## create a 10 x 10 grid of plot cells
par(mfrow=c(10,10))
## graph 100 stochastic map trees, sampled evenly from
## our set of 1,000
null<-sapply(mtrees[seq(1,100,by=1)],plot,colors=cols,lwd=1,ftype="off")

## compute posterior probabilities at nodes
pd<-summary(mtrees)
pd

## create a plot showing PP at all nodes of the tree
plot(pd,colors=cols,fsize=0.4,ftype="i",lwd=2,offset=0.4,ylim=c(-1,Ntip(anolis.tree)),cex=c(0.5,0.3))
## add a legend
legend("bottomleft",legend=levels(lnSVL),pch=22,pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)

## set margins
par(mar=c(5.1,4.1,2.1,2.1))
## graph marginal ancestral states and posterior
## probabilities from stochastic mapping
plot(fit.marginal$states[,1],pd$ace[1:anolis.tree$Nnode],pch=21,cex=1.2,bg="grey",xlab="Marginal scaled likelihoods",ylab="Posterior probabilities",bty="n",las=1,cex.axis=0.8)
lines(c(0,1),c(0,1),col="blue",lwd=2)

## create a "densityMap" object
anolis.densityMap<-densityMap(mtrees,states=levels(lnSVL)[2:1],plot=FALSE)

## update color gradient
anolis.densityMap<-setMap(anolis.densityMap,cols[2:1])
## plot it, adjusting the plotting parameters
plot(anolis.densityMap,fsize=c(0.3,0.7),lwd=c(3,4))
