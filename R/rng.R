#--functions for making and analysing RNG from 2d-binning spaces

prep.for.rng<-function(bDf,minl)
{
 #--take a 3-column (cov,len,gc) data.frame ('bDf'), subset by length ('minl') and return a matrix of log10(cov) and gc data for use in 'make.rng.from.2d.binning.space'
 bDfCut<-bDf[which(bDf$len>minl),]
 res<-cbind(log10(bDfCut$cov),bDfCut$gc)
 rownames(res)<-rownames(bDfCut)
 colnames(res)<-c("lgcov","gc")
 return(res)
}

make.rng.from.2d.binning.space<-function(bM)
{
 #--'bM' is a matrix with rows indexing contigs and 2 columns that index binning variables (e.g. two coverage vectors or one coverage and one GC, etc)
 bMColScale<-scale(bM,T,T)
 bMColScale.ppp<-spatstat::ppp(bMColScale[,1],bMColScale[,2],window=spatstat::owin(xrange=range(bMColScale[,1]),yrange=range(bMColScale[,2])))
 bMColScale.ppp.rng<-spatgraphs::spatgraph(bMColScale.ppp,type="RNG")
 bMColScale.ppp.rng.ig<-igraph::graph_from_adj_list(bMColScale.ppp.rng$edges)
 bMColScale.ppp.rng.ig.simple<-igraph::simplify(bMColScale.ppp.rng.ig)
 igraph::V(bMColScale.ppp.rng.ig.simple)$name<-rownames(bM)
 return(bMColScale.ppp.rng.ig.simple)
}

#make.rng.from.2d.binning.space.v1<-function(bM,...)
#{
# #--'bM' is a matrix with rows indexing contigs and 2 columns that index binning variables (e.g. two coverage vectors or one coverage and one GC, etc)
# bMColScale<-scale(bM,T,T)
# bMColScale.ppp.rng.ig<-rng(bMColScale,...)
# bMColScale.ppp.rng.ig.simple<-igraph::simplify(bMColScale.ppp.rng.ig)
# V(bMColScale.ppp.rng.ig.simple)$name<-rownames(bM)
# return(bMColScale.ppp.rng.ig.simple)
#}

#--extract k-neighbourhoods for all nodes...
extract.all.neighbourhoods<-function(rng,k)
{
 res<-igraph::ego(rng,k,igraph::V(rng),mindist=0)
 res<-lapply(res,FUN=function(x){x$name})
 names(res)<-igraph::V(rng)$name
 return(res)
}

#--summarise k-neighbourhoods...
summarise.neighbourhoods<-function(nhMembershipList,bDf)
{
 #--extract some properties of neighbourhoods and store in a table
 #--neighbourhood size
 nh.size<-sapply(nhMembershipList,length)
 
 #--and
 nh.mp=NULL
 nh.sd=NULL
 for(curnh in names(nhMembershipList))
 {
  curdata<-bDf[nhMembershipList[[curnh]],]
  nh.mp<-rbind(nh.mp,apply(curdata,2,mean))
  nh.sd<-rbind(nh.sd,apply(curdata,2,sd))
 }
 
 res<-data.frame(nh.size,nh.mp,nh.sd)
 colnames(res)<-c("size","lgcov.mean","gc.mean","lgcov.sd","gc.sd")
 rownames(res)<-names(nhMembershipList)
 
 return(res)
}
