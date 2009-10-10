####################### Internal Functions######################

######################################################
# The internal functions in LiquidAssociation Package
#######################################################

## General methods accessor
qqnorm2	<-function(object){
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect<-object
			n<-length(myvect)
			rmyvect<-rank(myvect)/(n+1)
			nmyvect<-sapply(rmyvect, qnorm)
			return(nmyvect)
}



stand<-function(object) {
			if (sum(is.na(object))>0)
			stop("Input vector contains missing value!!")
			myvect<-object
			ans<-(myvect-mean(myvect))/sd(myvect)
			return(ans)
}





form<-function(data){
		pts<-nrow(data)
		groupid<-rep(1:pts, each=2)
		y<-matrix(t(data[,1:2]),ncol=1, nrow=2*pts)
		i1<-rep(c(1,0), pts)
		i2<-rep(c(0,1), pts)
		x31<-rep(data[,3], each=2)*i1
		x32<-rep(data[,3], each=2)*i2
		visit<-rep(c(1,2), pts)
		dat<-as.data.frame(cbind(groupid, y, i1, i2,x31, x32, visit))
		colnames(dat)<-c("groupid", "y", "i1", "i2", "x31", "x32", "visit")
		return(dat)
}


##########################

boots.se<-function(object, num=1, dim){
		stand.data<-object
		n<-nrow(stand.data)
		sam<-sample(1:n, n, replace=TRUE)
		ndata<-stand.data[sam,]
		ndata[,dim]<-qqnorm2(ndata[,dim])
		nstand.ndata<-apply(ndata, 2, stand)
		ans<-LA(nstand.ndata)
		return(ans)
}


permute3.se<-function(object){
		stand.data<-object
		null<-sample(stand.data[,3])	
		ndata<-cbind(stand.data[,1:2], null)
		return(ndata)
}



boots.gla.se2<-function(object, num=1 ,cut=4, dim=3){
		stand.data<-object
		n<-nrow(stand.data)
		sam<-sample(1:n, n, replace=TRUE)
		ndata<-stand.data[sam,]
		ndata[,dim]<-qqnorm2(ndata[,dim])
		stand.ndata<-apply(ndata, 2, stand)
		ans<-GLA(stand.ndata, cut=cut, dim=dim)
		return(ans)
}
















