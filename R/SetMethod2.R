##############################
# Set S4 class functions
#
##############################


#### setMethods: GLA
setGeneric("GLA", function(object, ...) standardGeneric("GLA"))
setMethod("GLA", signature(object="eSet"), function(object, geneMap=NULL, ...){
		dat<-t(exprs(object))
		GLA(dat, geneMap,...)
	}
)
setMethod("GLA", signature(object="matrix"), 
	function(object, cut=4, dim=3, geneMap=NULL){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		data<-object[!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3]),]
		data[,dim]<-qqnorm2(data[,dim])
		data<-apply(data, 2, stand)
		x<-seq(0,1, length=cut)
		br<-quantile(data[,dim], prob=x)
		index<-as.numeric(cut(data[,dim], breaks=br))
		tab<-table(index)
		tab<-tab[tab > 2]
		vect<-as.numeric(names(tab))
		m2<-rep(0, length(vect))
		cor.e<-rep(0, length(vect))
		gla.x2<-rep(0, length(vect))
		for ( i in 1:length(vect)){
			p<-which(index==vect[i])
			m2[i]<-mean(data[p,dim])
			cor.e[i]<-cor(data[p, -dim])[1,2]
			gla.x2[i]<-cor.e[i]*m2[i]
		}
		ans<-mean(gla.x2)
		names(ans)<-paste("GLA(", colnames(object)[1], ",", colnames(object)[2], "|", colnames(object)[3],")", sep="")
		return(ans)
	 }
)

		
		
###############

setGeneric("LA", function(object,...) standardGeneric("LA"))
setMethod("LA", signature(object="eSet"), function(object, geneMap=NULL,...){
		dat<-t(exprs(object))
		LA(dat, geneMap,...)
	}
)
setMethod("LA", signature(object="matrix"), 
	function(object, geneMap=NULL,dim=3){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		ndat<-object[!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3]),]
		ndat[,dim]<-qqnorm2(ndat[,dim])
		ndat<-apply(ndat,2, stand)
		anss<-apply(ndat,1,prod)
 		ans<-mean(anss)
 		names(ans)<-paste("LA(", colnames(object)[1], ",", colnames(object)[2], ",", colnames(object)[3], ")", sep="")
		return(ans)
	}
)


setClass("CNM", representation(Model="character", output="matrix"))
setMethod("print","CNM",
          function(x){
          	cat("\n")
          	cat(x@Model);cat("\n")
          	print(x@output)
          	}
)
setMethod("show", signature(object="CNM"),
          function(object){
          	cat(object@Model);cat("\n")
          	print(object@output)
          	}
)

setGeneric("CNM.full", function(object,...) standardGeneric("CNM.full"))
setMethod("CNM.full", signature(object="eSet"), function(object, geneMap=NULL, ...){
		dat<-t(exprs(object))
		CNM.full(dat, geneMap, ...)
	}
)
setMethod("CNM.full", signature(object="matrix"), 
	function(object, geneMap=NULL, dim=3){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		dat<-object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
		dat[,dim]<-qqnorm2(dat[,dim])
		dat<-apply(dat, 2, stand)
	    data<-form(dat)
		zcor<-cbind(1,dat[,3])
		fit<-geese(y ~  x31 + x32 -1, id=groupid, data=data, sformula= ~ i1 + i2 +  x31 + x32 -1, waves=visit, 
		sca.link="log", cor.link="fisherz", corstr="userdefined", control=geese.control(trace=FALSE),
		zcor=zcor, jack=TRUE, j1s=TRUE, fij=TRUE)
		out1<-summary(fit)
		tabout<-matrix(0, nrow=8, ncol=4)
		tabout[,1]<-c(out1[[3]]$estimate[1], out1[[3]]$estimate[2], out1[[2]]$estimate[1], out1[[1]]$estimate[1], out1[[1]]$estimate[2], out1[[3]]$estimate[3], out1[[3]]$estimate[4],  out1[[2]]$estimate[2])  
		tabout[,2]<-c(out1[[3]]$san.se[1], out1[[3]]$san.se[2],out1[[2]]$san.se[1], out1[[1]]$san.se[1], out1[[1]]$san.se[2],out1[[3]]$san.se[3], out1[[3]]$san.se[4], out1[[2]]$san.se[2])
		tabout[,3]<-c(out1[[3]]$wald[1], out1[[3]]$wald[2],out1[[2]]$wald[1], out1[[1]]$wald[1], out1[[1]]$wald[2],out1[[3]]$wald[3], out1[[3]]$wald[4], out1[[2]]$wald[2])
     	tabout[,4]<-c(out1[[3]]$p[1], out1[[3]]$p[2],out1[[2]]$p[1], out1[[1]]$p[1], out1[[1]]$p[2],out1[[3]]$p[3], out1[[3]]$p[4],  out1[[2]]$p[2])
		colnames(tabout)<-c("estimates", "san.se", "wald", "p value")
		rownames(tabout)<-c("a3", "a4", "a5", "b1", "b2", "b3", "b4", "b5")
		Model<-paste("Model: CNM(", colnames(object)[1], ",", colnames(object)[2], "|", colnames(object)[3],")", sep="")
	new("CNM", Model=Model, output=tabout)
	}
)


setGeneric("CNM.simple", function(object,...) standardGeneric("CNM.simple"))
setMethod("CNM.simple", signature(object="eSet"), function(object, geneMap=NULL,...){
		dat<-t(exprs(object))
		CNM.simple(dat, geneMap,...)
	}
)
setMethod("CNM.simple", signature(object="matrix"), 
	function(object, geneMap=NULL,dim=3){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		dat<-object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
		dat[,dim]<-qqnorm2(dat[,dim])
		dat<-apply(dat,2,stand)
		data<-form(dat)
		zcor<-cbind(1,dat[,3])
		geefit<-geese(y ~ 1, id=groupid, data=data, waves=visit, 
		sca.link="log", cor.link="fisherz", corstr="userdefined", 
		zcor=zcor, jack=TRUE, j1s=TRUE, fij=TRUE, control=geese.control(trace=FALSE))
		out<-summary(geefit)
		tabout<-matrix(0, nrow=4, ncol=4)
		tabout[,1]<-c(out[[3]]$estimate[1], out[[3]]$estimate[1], out[[2]]$estimate[1], out[[2]]$estimate[2])
		tabout[,2]<-c(out[[3]]$san.se[1],   out[[3]]$san.se[1],   out[[2]]$san.se[1],   out[[2]]$san.se[2])
		tabout[,3]<-c(out[[3]]$wald[1],     out[[3]]$wald[1],     out[[2]]$wald[1],     out[[2]]$wald[2])
    	tabout[,4]<-c(out[[3]]$p[1],        out[[3]]$p[1],       out[[2]]$p[1],        out[[2]]$p[2])
		colnames(tabout)<-c("estimates", "san.se", "wald", "p value")
		rownames(tabout)<-c("a3", "a4", "a5", "b5") 
		 Model<-paste("Model: CNM(", colnames(object)[1], ",", colnames(object)[2], "|", colnames(object)[3],")", sep="")
	new("CNM", Model=Model, output=tabout)
	}
)





####### Hypothesis Testing methods accessor

setGeneric("getsGLA", function(object, ...) standardGeneric("getsGLA"))
setMethod("getsGLA", signature(object="eSet"), function(object, geneMap=NULL,...){
		dat<-t(exprs(object))
		getsGLA(dat, geneMap,...)
	}
)
setMethod("getsGLA", signature(object="matrix"),
	function(object, boots=30, perm=100, cut=4, dim=3, geneMap=NULL){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		dat<-object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
		ans<-matrix(0, nrow=perm, ncol=2)
		dat[,dim]<-qqnorm2(dat[,dim])
		nsdata<-apply(dat, 2, stand)
		gla1<-GLA(nsdata, cut=cut, dim=dim)
		bootsse1<-sapply(1:boots, boots.gla.se2, object=nsdata, cut=cut, dim=dim)
 		bse1<-sd(bootsse1)
		test1<-gla1/bse1
		for ( i in 1:perm){
			datanull<-permute3.se(nsdata)
			gla2<-GLA(datanull, cut=cut, dim=dim)
			bootsse2<-sapply(1:boots, boots.gla.se2, object=datanull, cut=cut, dim=dim)
 	 		bse2<-sd(bootsse2)
			ans[i]<-gla2/bse2
		}
		pvalue<-2*sum(abs(test1) <= ans)/perm
		out<-c(test1, pvalue)
		names(out)<-c("sGLA", "p value")
		return(out)
		}
)



setGeneric("getsLA", function(object, ...) standardGeneric("getsLA"))
setMethod("getsLA", signature(object="eSet"), function(object, geneMap=NULL,...){
		dat<-t(exprs(object))
		getsLA(dat, geneMap,...)
	}
)
setMethod("getsLA", signature(object="matrix"),
	function(object, boots=30, perm=100, geneMap=NULL, dim=3){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		ans<-rep(0, perm)
		data<-object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
     	data[,dim]<-qqnorm2(data[,dim])
     	sdata<-apply(data, 2, stand)
      	la1<-LA(sdata)
      	bootsse1<-sapply(1:boots, boots.se, object=sdata)
      	bse1<-sd(bootsse1)
      	test1<-la1/bse1
		for ( i in 1:perm){
			datanull<-permute3.se(sdata)
			la2<-LA(datanull)
			bootsse2<-sapply(1:boots, boots.se, object=datanull, dim=dim)
        	bse2<-sd(bootsse2)
			ans[i]<-la2/bse2
		}
		pvalue<-2*sum(abs(test1) <= ans)/perm
		out<-c(test1, pvalue)
		names(out)<-c("sLA", "p value")
		return(out)
	}
)


###### plot function

setGeneric("plotGLA", function(object, ...) standardGeneric("plotGLA"))
setMethod("plotGLA", signature(object="eSet"), function(object, geneMap=NULL,...){
		dat<-t(exprs(object))
		plotGLA(dat, geneMap,...)
	}
)
setMethod("plotGLA", signature(object="matrix"),
	function(object, cut=4, dim=3, filen="b1b2dat5", save=FALSE, geneMap=NULL,...){
		if (!is.numeric(object))
			stop("Input matrix must be numeric")
		if(ncol(object)!=3)
			stop("Input data must have three variables")
		 if (length(colnames(object))!=3 & length(geneMap)!=3)
			stop("Please specify the names of three variables")
		if (length(geneMap)==3)
			colnames(object)<-names(geneMap)
		data<-object[(!is.na(object[,1]) & !is.na(object[,2]) & !is.na(object[,3])),]
		dimm<-c(1,2,3)[-dim]
		x<-seq(0,1, length=cut)
		br<-quantile(data[,dim], prob=x)
		index<-as.numeric(cut(data[,dim], breaks=br, include.lowest=TRUE))
		tab<-table(index)
		tab<-tab[tab > 2]
		vect<-as.numeric(names(tab))
		x3name<-colnames(data)[dim]
		for ( i in 1:length(vect)){
			p<-which(index==vect[i])
			filename<-paste(filen, dim, i, ".pdf", sep="")
			corr<-cor(data[p, dimm[1]], data[p, dimm[2]])
			minx<-min(data[p, dim])
			maxx<-max(data[p, dim])
			lx<-min(data[p, dimm[1]])
			ly<-max(data[p, dimm[2]])
			if (!save){
				plot(data[p, dimm], main=paste(round(minx,1), "<", x3name,"<", round(maxx,1), sep=""), ...)
				pt<-colMeans(data[p, dimm])
				points(pt[1],pt[2], pch=16, cex=2, col=3)
				abline(h=0,lty=2, col="lightgrey")
				abline(v=0, lty=2, col="lightgrey")
			legend(lx+0.2, ly-0.2, paste("corr=", round(corr,2)), text.col=2,cex=2, bty="n")
				Sys.sleep(0.5)
			}else{
				pdf(filename)
				plot(data[p, dimm], main=paste(round(minx,1), "<X", dim,"<", round(maxx,1), sep=""), ...)
				pt<-colMeans(data[p, dimm])
				points(pt[1],pt[2], pch=16, cex=2, col="grey")
				abline(h=0,lty=2, col="lightgrey")
				abline(v=0, lty=2, col="lightgrey")
			legend(lx+0.2, ly-0.2, paste("corr=", round(corr,2)), text.col=2, cex=2, bty="n")
				dev.off()
			}
		}
	}
)

