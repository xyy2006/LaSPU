###################### 5 tests for association, possibly including interactions
######################   in the possible presence of nuisance parameters:
######################        Score, SSU, SSUw, Sum, UminP tests.
###########References:  Pan W, Basu S, Shen X (2011). Adaptive Tests for Detecting Gene-Gene and Gene-Environment Interactions.

######################Wei Pan, Division of Biostatistics, U of Minnesota
######################Email: weip@biostat.umn.edu
#### V1, Feb 16, 2011


library("mvtnorm")

#Input:
#      Y: a vector of 0-1 response, nx1;
#      X: design matrix for covariates (e.g. gene-gene interactions) of
#           interest (i.e. to be tested), nxk;
#      Z: design matrix for nuisance covariates, nxk2;
#Goal:
#      fit model: logit(Pr(Y=1)) = a0+ a1*X+a2*Z ,
#      test H0: a1=0.
#Output: p-values for Score, SSU, SSUw, Sum, UminP

#Examples: suppose that one has two groups of variables in design matrices
#          X1 and X2 respectively, and Y is a binary response vector,
#          1) if one would like to test the marginal effects of X1, then
#             Tests5V1(Y, X=X1);
#          2) if one would like to test the effects of X1 after ajusting for X2,
#             Tests5V1(Y, X=X1, Z=X2);
#          3) if one would like to test the effects of X1 in the presence of 
#             possible interactions with X2, then
#             Tests5V1(Y, X=cbind(X1, X12), Z=X2), 
#             where X12 is the design matrix for interactions;
#          4) if one would like to test on the interaction between X1 and X2
#             after adjusting for the main effects of X1 and X2, then
#             Tests5V1(Y, X=X12, Z=cbind(X1,X2)), 
#             where X12 is the design matrix for interactions;

############################################
############################################
############################################

Tests5V1<-function(Y, X, Z=NULL ){

n<-length(Y)
k<-ncol(X)
k2<-ncol(Z)

#######construction of the score vector and its cov matrix:
if (is.null(Z)){
## NO nuisance parameters:
   Xg <- X
   Xbar<-apply(Xg, 2, mean)
   Xgb<-Xg
   for(i in 1:nrow(Xg))
      Xgb[i,]<-Xg[i,]-Xbar
   ##score vector:
   U<-t(Xg) %*% (Y-mean(Y))
   ##cov of the score stats:
   CovS<- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)

} else {
## with nuisance parameters:
tdat1<-data.frame(trait=Y, Z)
fit1<-glm(trait~.,family="binomial",data=tdat1)
pis<-fitted.values(fit1)
Us<-matrix(0, nrow=n, ncol=k)
for(i in 1:k){
    tdat2<-data.frame(X1=X[,i], Z)
    fit2<-glm(X1~.,data=tdat2)
    X1mus<-fitted.values(fit2)
    Us[, i]<-(Y - pis)*(X[,i] - X1mus)
    }
U<-apply(Us, 2, sum)
CovS<-matrix(0, nrow=k, ncol=k)
for(i in 1:n)
  CovS<-CovS + Us[i,] %*% t(Us[i,])
}

pSco<-Score(U, CovS)
pSSU<-SumSqU(U, CovS)
pSSUw<-SumSqUw(U, CovS)
pSum<-Sum(U, CovS)
pUminP<-UminP(U, CovS)

return(c(pSco, pSSU, pSSUw, pSum, pUminP))

}  

SumSqU<-function(U, CovS){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTg1<-as.numeric(1-pchisq(Tscore, 1))
        }
     else {
     #it's possible U=0 and Cov(U)=0:
     if (all(abs(U)<1e-20)) pTg1<-1 else{
     Tg1<- t(U) %*% U
     ##distr of Tg1 is sum of cr Chisq_1:
       cr<-eigen(CovS, only.values=TRUE)$values
     ##approximate the distri by alpha Chisq_d + beta:
       alpha1<-sum(cr*cr*cr)/sum(cr*cr)
     beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
     d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
     # alpha1<-as.real(alpha1)
     # beta1<-as.real(beta1)
     # d1<-as.real(d1)
     alpha1<-as.numeric(alpha1)
     beta1<-as.numeric(beta1)
     d1<-as.numeric(d1)
     pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
     }
     }
return(pTg1)
}

SumSqUw<-function(U, CovS){
     ###STRANGELY: diag(A) fails if A is 1X1:
# diag(CovS[40,40])
#<0 x 0 matrix>
#> CovS0<-CovS[40,40]
#> CovS0
#[1] 0.24975
#> diag(CovS0)
#<0 x 0 matrix>
#> diag(1)
#     [,1]
#[1,]    1

     if (is.null(dim(CovS))) {# only one-dim:
        Tg2<- sum(U^2 /CovS)
        if (is.na(Tg2) || is.infinite(Tg2) || is.nan(Tg2)) Tg2<-0
        pTg2<-as.numeric(1-pchisq(Tg2, 1))
        }
     else {
     #it's possible U=0 and Cov(U)=0:
     if (all(abs(U)<1e-20)) pTg2<-1 else{
     diagCovS<-diag(CovS)
     diagCovS<-ifelse(diagCovS>1e-10, diagCovS, 1e-10)
     #Tg2<- t(U) %*%  diag(1/diagCovS) %*% U
     Tg2<- sum(U^2 /diagCovS)
     ##distr of Tg1 is sum of cr Chisq_1:
       cr<-eigen(CovS %*% diag(1/diagCovS), only.values=TRUE)$values
     ##approximate the distri by alpha Chisq_d + beta:
       alpha2<-sum(cr*cr*cr)/sum(cr*cr)
       beta2<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
       d2<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
       alpha2<-as.numeric(alpha2)
       beta2<-as.numeric(beta2)
       d2<-as.numeric(d2)
     pTg2<-as.numeric(1-pchisq((Tg2-beta2)/alpha2, d2))
     }
     }
   pTg2
}

 ##########score test:
Score<-function(U, CovS){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTscore<-as.numeric(1-pchisq(Tscore, 1))
        }
    else {
    ##gInv of CovS:
      CovS.edecomp<-eigen(CovS)
      CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
      inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
      P<-solve(CovS.edecomp$vectors)
      gInv.CovS<-t(P) %*% diag(inveigen) %*% P
    Tscore<- t(U) %*% gInv.CovS  %*% U
    pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )
    }
pTscore
}

##########SumTest########################
Sum<-function(U, CovS){
     #it's possible U=0 and Cov(U)=0:
     if (all(abs(sum(U))<1e-20)) pTsum<-1 else{
      a<-rep(1, length(U))
      Tsum<- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
      pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
      }
pTsum
}

##########UminP Test########################
UminP<-function(U, CovS){

     if (is.null(dim(CovS))) {# only one-dim:
        Tu<- sum(U^2 /CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) Tu<-0
        pTu<-as.numeric(1-pchisq(Tu, 1))
        }
     else{
          ####it's POSSIBLR Ui=0 and CovS[i,i]=0!
          Tu<-as.vector(abs(U)/(sqrt(diag(CovS)) + 1e-20) )
    k<-length(U)
    V <- matrix(0,nrow=k, ncol=k)
    for(i in 1:k){
      for(j in 1:k){
          if (abs(CovS[i,j])>1e-20)
             V[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
          else   V[i,j] <- 1e-20
        }
      }
    pTu <- as.numeric(PowerUniv(Tu,V))	# call MVN to get the p-value, may not be good for RVs.
        }

pTu
}

#-------------------------------#
# a subfunc using U vector and V - Vcov matrix, according to mvnorm dist, return the pvalue.
PowerUniv <- function(U,V){
      n <- dim(V)[1]

      x <- as.numeric(max(abs(U)))
      TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))

      return(TER)
}
