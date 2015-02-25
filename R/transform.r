#################################################
## These are the functions to transform the matrix to the diagonalized design matrix for GEE
## created by Yiwei Zhang
### tranform.x is to transform x
## INPUT:
# dat is n by p matrix, with n observation and p variables
# gap is the number of traits per observation (q)
## OUTPUT:
# a n*gap by p*gap design matrix
### transform.y is to transform y
##INPUT: y is a n by q matrix for n observation and q  traits (which is gap)
##OUTPUT: a vector of n*q



transform.x<-function(dat, gap){
  dat=as.matrix(dat)
  n.col=ncol(dat)*gap
  n.row=nrow(dat)*gap
  k=ncol(dat)

dat.new=matrix(0,ncol=n.col, nrow=n.row)
for (r in 1:gap) dat.new[seq(r,n.row,gap),(k*(r-1)+1):(k*r)]=as.matrix(dat)

return(dat.new)
}

transform.y<-function(y,gap){
 if (!is.null(dim(y))) n=nrow(y) else n=length(y) 
newy=matrix(0,nrow=n*gap,ncol=1)
for (r in 1:gap) newy[seq(r,n*gap,gap),]=y[,r]
return(newy)
}


