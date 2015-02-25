###### expedited version of gee.spu.r, created by Yiwei Zhang, updated by Yang Yang.
###### It is to calculate p-value of aSPU, aSPU.w and aSPU.sco with permutation based P-value (mainly fo rare variant scenario) between subject but keep
###### the time order for longitudinal measurements. Score test within updated; provide a Paspu_aspuw.sco, which combines aSPU, aSPUw and score statistics by minimum-pvalue combining method.
###### Because when B=10to9, R can not store the whole matrix B by length(gamma)+1
###### So here we draw same samples for each gamma and calculate p-value for gamma=1,2,3,...,Inf
###### To calculate aSPU, for gamma(j), the distribution of the samples are calculated (minp0=P0s(j)),and for gamma(j+1), if P0s(j+1) is smaller than P0s(j), than the values are replaced (minp0=P0s(j+1))
###### aSPU p-valu is calculated as proportion of min(P)>minp0
###### wSPU added on 9/7/2013
###### paralleled on 01/28/2014
###################################################################################
# usage: source this script in main script.
# Note: foreach backend should be registered in main script, e.g. doMC, doSNOW, etc.
#		for this script, suggest always enable .parallel_overGamma, it saves time (since we have 9 lambda value to parallel)
#		2, we put Uminp caculated with gamma=0~Inf, since only need one time, we put it with 1st gamma encoutered.
#		3, suggest .parallel_overGamma and .parallel_overB cannot be both TRUE (avoid oversubscrib of parallel resource.)
#		4, if use permutation, suggest .parallel_overP  = TRUE; other wise, very slow!!!
###################################################################################
# argument definitions:
# B: the times for either simulation or permutation.
#======================================#
# original updated.
SPU3_withScoGee_RV_version<-function(U,V,gamma=c(0:8,Inf),B, .seed=12345,.parallel_overGamma=FALSE, .parallel_overB = FALSE, cores_overB = ifelse(.parallel_overB, 2, 1), usePermutationNullU = FALSE, permutationNullU){
	# if (!.parallel) registerDoSEQ()	# disable parallel backend registered at parent.frame(), which called this function.
	if (!.parallel_overGamma) `%dopar%` <- `%do%`	# change local definition within function.
	
	a=ifelse(diag(V)>1e-16,diag(V),1e-16)	# for digit accuracy consideration, to avoid 1/0 case.
	# a=ifelse((V)>1e-16,(V),1e-16)	# use full V instead of diagnol V.
	# a=diag(V)	# a fast way to equate a to diagnol of V, however, if 0 exists in V, NA will return.
	##observed spu
	spu=spuw=rep(NA,length(gamma))
	U_div_by_sqrta <- U/sqrt(a)
	#-------calculate the spu and spuw statistics--------------------#
	for (g in 1:length(gamma)){
		if (gamma[g]==0) {
			# spu[g]<- spuw[g] <- t(U)%*%ginv(V)%*%U # put score test statistics at the 1st index place.
			spu[g]<- spuw[g] <- Score_robust(U, CovS=V)$Score.statistics # put score test statistics at the 1st index place.
		}
		if ( gamma[g]>0 & gamma[g]<Inf ) {
			spu[g]<-sum(U^gamma[g])
			spuw[g]<-sum(U_div_by_sqrta^gamma[g])
		} 
		if (gamma[g]==Inf) {
			spu[g]<-max(abs(U))
			spuw[g]<-max(abs(U_div_by_sqrta))
		} 
	}	# now we will have a g here.
	#-------------------------------#
	T<-max(U^2/a)	# UminP observed T statistics
	cat("statistic calculated","\n")
	## perm 
	p=pw=rep(NA,length(gamma))
	# T1s=T1sw=numeric(B)
	# count=0
	# gamma_statistics_return <- list()
	minp0=minp0w=0
  
	# str(Results)
	# clusterSetupRNGstream(cl,seed=rep(seed,6))
	# clusterSetupRNG(cl)
	#--------------------------#
	# the func to do get simulation based statistics.
	gamma_statistics_return <- function (b=1,g) {
		u.null<-uNull[b,]	# extract a row here, uNull is in calling env.
		u.null_div_by_sqrta <- u.null/sqrt(a)
		if (gamma[g]==0) {
			# T1s<- T1sw <- t(u.null)%*%ginv(V)%*%u.null # put score test statistics at the 1st index place, slower than g=other values, i.e. gamma != 0, due to matrix multiplication operations.
			T1s<- T1sw <- Score_robust(u.null, CovS=V)$Score.statistics # put score test statistics at the 1st index place, slower than g=other values, i.e. gamma != 0, due to matrix multiplication operations.
			#---not better at all for above---#
			# T1s<- T1sw <- crossprod(u.null,ginv(V) ) %*% u.null # put score test statistics at the 1st index place, slower than g=other values, i.e. gamma != 0.
		}
		if ( gamma[g]>0 & gamma[g]<Inf ) {
			T1s<-sum(u.null^gamma[g])
			T1sw<-sum(u.null_div_by_sqrta^gamma[g])
		} 
		if (gamma[g]==Inf) {
			T1s<-max(abs(u.null))
			T1sw<-max(abs(u.null_div_by_sqrta))
		} 
		  
		if (gamma[g]==1){	# random g, just do it once for uminP in this function. not on every g. uminP not related to any g values.
		  t0<-max(u.null^2/a); # to be later compared with T in #.41
		  # count <- as.integerI(t0>T)
		} else t0 <- NA	# for uniform umin/minP.
	   return(c(b = b, T1s = T1s, T1sw = T1sw,t0 = t0))
	}
	# and need the matching uNull.
	
	if (!usePermutationNullU) {
		set.seed(.seed) # to ensure the same samples are drawn for each gamma. i.e. set same seed on each slave. To set a seed on a series of slaves please use clusterSetupRNG. So batch result will be the same but not within each slave.
		uNull <- mvrnorm(B,rep(0,length(U)),V)	# simulate uNull should come from multivariate N(0,V), only right for common variants.
	} else uNull <- permutationNullU	# need to be in the same # of row as mvrnorm(B,...)
	# OR use permutation based
	#--------------------------------------------#

  
  # start computation for each gamma.
	# Results <- foreach (g = rep(1,9),.packages="MASS",.verbose=FALSE) %dopar% {	# # for debug use all same g.
	#######################################################
	##############start computation for each gamma#########
	#######################################################
	#------------------now get the aSPU(w).sco by using all gammas.-----------------------------------#
	Results <- foreach (g = 1:length(gamma),.packages=c("MASS","doSNOW","parallel"),.verbose=F,.errorhandling="pass", .export = c('Score_robust')) %dopar% {	# 'Score_robust' is not defined in local environment.
  # for (g in 1:length(gamma)) {	# now this loop level are serial for purpose.
		# collect_garbage()
		# catt(g,'started!')
		# set.seed(seed) # to ensure the same samples are drawn for each gamma. i.e. set same seed on each slave. To set a seed on a series of slaves please use clusterSetupRNG. So batch result will be the same but not within each slave.
		#
		# if (gamma[g] != 0) {
		if(!.parallel_overB) {
		  T_matrix_temp <-simplify2array( lapply(1:B,gamma_statistics_return, g=g) )	# the g take effect within gamma_statistics_return
		} else {
			# # or use fork parallel.
		  T_matrix_temp <-simplify2array( mclapply(1:B,gamma_statistics_return,g=g,mc.cores=cores_overB ) )	# the g take effect within gamma_statistics_return; 3 cores * 12 = 36 working on a 12 cores node.
		}
		#
		# names(gamma_statistics_return) <- gamma
		################################################################
		#-------calculating the pvalue based on B+1 as denominator-----#
		# "adding 1" in the numerator essentially treats the
		# observed test stat based on the original data as one
		# from the null distribution under H0. again it is
		# clear that this "adding 1" does not make any difference
		# as B goes to infinity.
		p= mean( abs(spu[g]) <= c( abs(spu[g]), abs(T_matrix_temp["T1s",])) )	# the spu(g) pvalue for a specific gamma. the observed one.
		# obs <= simulated, mean we simulate the sample more extreme or as extreme as the observed one.
		# Note that there should be an equality sign in
		# "sum I(|Tb| >= |Tobs| + 1 )/B",
		# NOT "sum I(|Tb| > |Tobs| + 1 )/B",
		# which is important for discrete traits and RVs.

		pw= mean( abs(spuw[g]) <= c(abs(spuw[g]),abs(T_matrix_temp["T1sw",])) ) # the observed one pvalue for SPUw(g).
		
		#
		if (g==1) t0_vec <- T_matrix_temp["t0",]  else t0_vec <- NULL	# for UminP statistics.
		
		#--------calculating the pvalue for each b in B (one sample ) within the same power gamma.-----------#
		# a<- 1:1000
		# a[998] = a[999]	# so there is a tie.
		# mean(a[-999]>= a[999])
		# [1] 0.002002002
		# mean(a[-998]>= a[998])
		# [1] 0.002002002
		# (1000 - rank(abs(a),ties.method="min"))/ (1000 - 1)
		# 0.002002002 0.002002002 0.000000000
		# and then if we count B-1+1.
		# > mean(a[]>= a[999])
		# [1] 0.003
		# > mean(a[]>= a[998])
		# [1] 0.003
		# > mean(a[]>= a[997])
		# [1] 0.004
		# > mean(a[]>= a[999])
		# [1] 0.003
		# > mean(a[]>= a[1000])
		# [1] 0.001
		# same as
		# (1000 - rank(abs(a),ties.method="min") + 1 )/ (1000)
		P0s=(B-rank(abs(T_matrix_temp["T1s",]),ties.method="min") + 1) / (B)	# consider B-1+1
		P0sw=(B-rank(abs(T_matrix_temp["T1sw",]),ties.method="min") + 1) / (B)	# consider B-1+1
		# validated by:
		# mean(a[-998]>= a[998]) = 0.002002002
		
		# return(list(minp0 = minp0, minp0w = minp0w, p = p, pw = pw ))
		# catt(g,'finished!')
		return(list(t0_vec = t0_vec,P0s = P0s, P0sw = P0sw, p = p, pw = pw ))
	}	# for each g
	names(Results) <- gamma	# this obj may not need.
	# browser()
  #	#######################################################
  #	#######################################################
  #	#######################################################
	cat("all g finished!\n")
	# browser()
  	#-------------get the minimal SPU(gamma) across every b--------#
	# foreach( g = 1:length(Results)) %do% {	# avoid parallel here, too easy task will cost more time on communication then computation.
		# temp <- Results[[g]]
		# if (g==1) minp0=temp$P0s else minp0[which(minp0>temp$P0s)] <- temp$P0s[which(minp0>temp$P0s)]
		# if (g==1) minp0w=temp$P0sw else minp0w[which(minp0w>temp$P0sw)] <- temp$P0sw[which(minp0w>temp$P0sw)]
		# cat("g=",g,"finished!\n")
	# }
	# equivalent to below:
	Results_excludeScoGeeIndex <- Results[- match(0,names(Results) ) ] # need to guarantee the gamma = 0 is used to represent the geescore Test statistics.
	#----------------------------------#
	### 1, aSPU and aSPU.w with geescore test statistics involved,
	### 	so produce aSPU(.w).sco
	minp0.sco <- apply( sapply(Results, function(x) x$P0s),1, min) # find for each b in B, that is a sample, the minimum pvalue across different power gammas. p-aSPU(b), for b = 1,...,B
	minp0w.sco <- apply( sapply(Results, function(x) x$P0sw),1, min) # find for each b in B, that is a sample, the minimum pvalue across different power gammas, that is the p-aSPUw(b), for b = 1,...,B
	
	p.sco = sapply(Results, function(x) x$p)	# observed SPU pvalue across gammas
	pw.sco = sapply(Results, function(x) x$pw)
	Paspu.sco <- mean( c(min(p.sco),minp0.sco) <= min(p.sco) )	# consider B+1, count how many simulated is even smaller than obs.
	Paspuw.sco <- mean( c(min(pw.sco),minp0w.sco) <= min(pw.sco) )
	#---#
	# generate the vector containing separate p-values for SPU(g)s and aSPU. For function return usage.
	p_vector_output.sco <- c(p.sco,aspu.sco = Paspu.sco)
	pw_vector_output.sco <- c(pw.sco,aspu.w.sco = Paspuw.sco)
	#---#
	# the overall P-aSPU+aSPUw_score. For function return usage.
	Paspu_aspuw.sco <- mean( c( min(p.sco,pw.sco), pmin.int(minp0.sco,minp0w.sco) ) <= min(p.sco,pw.sco) )	# still B+1, but parallel min of the minp0.sco and minp0w.sco.
	
	#----------------------------------#
	### 2, aSPU and aSPU.w, i.e., excluding the score gee test statistics.
	minp0 <- apply( sapply(Results_excludeScoGeeIndex, function(x) x$P0s),1, min) # find for each b in B, that is a sample, the minimum pvalue across different power gammas.
	minp0w <- apply( sapply(Results_excludeScoGeeIndex, function(x) x$P0sw),1, min) # find for each b in B, that is a sample, the minimum pvalue across different power gammas.
	
	p = sapply(Results_excludeScoGeeIndex, function(x) x$p)	# SPU pvalue across gammas
	pw = sapply(Results_excludeScoGeeIndex, function(x) x$pw)
	Paspu <- mean( c(min(p),minp0) <= min(p) )	# consider B+1
	Paspuw <- mean( c(min(pw),minp0w) <= min(pw) )
	#
	p_vector_output <- c(p,aspu = Paspu)
	pw_vector_output <- c(pw,aspu.w = Paspuw)
	#---#
	# the overall P-aSPU+aSPUw. For function return usage.
	Paspu_aspuw <- mean( c( min(p,pw), pmin.int(minp0,minp0w) ) <= min(p,pw) )	# still B+1, but parallel min of the minp0.sco and minp0w.sco.
	
	
	#--------------------------------------#
	# umin= mean( c(Results[[1]][["t0_vec"]],T) >= T )	# consider B+1
	# or below one dont use g index computation about t0_vec any more. Preferred.
	t0 <- sapply(1:B, function(b) max(uNull[b,]^2/a ))
	uminP = mean( c(t0,T) >= T )	# consider B+1, T is observed uminP.
	
	
	
	return(list(spu=p_vector_output,spu.w=pw_vector_output,minP=uminP,spu.sco = p_vector_output.sco, spu.w.sco = pw_vector_output.sco, Paspu_aspuw.sco = Paspu_aspuw.sco, Paspu_aspuw = Paspu_aspuw))
}

######################################################################
Score_robust<-function(U, CovS){
	if (is.null(dim(CovS))) {# only one-dim:
		Tscore<- sum(U^2 /CovS)
		if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
		pTscore<-as.numeric(1-pchisq(Tscore, 1))
		return( data.frame(df=1, Score.statistics=Tscore, pval=pTscore ) )
	} else {
	##gInv of CovS:
	  CovS.edecomp<-eigen(CovS)
	  CovS.rank<-sum(abs(CovS.edecomp$values)> 1e-8)
	  inveigen<-ifelse(abs(CovS.edecomp$values) >1e-8, 1/CovS.edecomp$values, 0)
	  P<-solve(CovS.edecomp$vectors)
	  gInv.CovS<-t(P) %*% diag(inveigen) %*% P
	  Tscore<- t(U) %*% gInv.CovS  %*% U
	  pTscore<-as.numeric( 1-pchisq(Tscore, CovS.rank) )
	  return( data.frame(df=CovS.rank, Score.statistics=Tscore, pval=pTscore ) )
	}	
}




#########################################
# #=====================================#
# # with full V involved, repalce '/sqrt(a)' with '%*% with(svd(ginv(a) ),u %*% diag(sqrt(d)) )'. For pd symmetric matrix, ususally the case Rw is, with(svd(ginv(a) ),u %*% diag(sqrt(d)) ) %*% t(itself) = ginv(a).
# SPU3<-function(U,V,gamma=c(1:8,Inf),B, seed=12345,.parallel=FALSE){
	# if (!.parallel) registerDoSEQ()	# disable parallel backend registered at parent.frame(), which called this function.
	
	# # a=ifelse(diag(V)>1e-16,diag(V),1e-16)	# for digit accuracy consideration, to avoid 1/0 case.
	# a=ifelse((V)>1e-16,(V),1e-16)	# use full V instead of diagnol V.
	# # a=diag(V)	# a fast way to equate a to diagnol of V, however, if 0 exists in V, NA will return.
	# ##observed spu
	# spu=spuw=rep(NA,length(gamma))
	
	# #-------calculate the spu and spuw statistics--------------------#
	# for (g in 1:length(gamma)){
	# if (gamma[g]<Inf) spu[g]<-sum(U^gamma[g]) else spu[g]<-max(abs(U))
	# if (gamma[g]<Inf) spuw[g]<-sum((U %*% with(svd(ginv(a) ),u %*% diag(sqrt(d)) ) )^gamma[g]) else spuw[g]<-max(abs(U%*% with(svd(ginv(a) ),u %*% diag(sqrt(d)) )))
	# }
	
	# #------UminP observed T statistics-------------------------#
	# T<-max(U^2/a)
	# cat("statistic calculated","\n")
	# ## perm 
	# p=pw=rep(NA,length(gamma))
	# # T1s=T1sw=numeric(B)
	# # count=0
	# # gamma_statistics_return <- list()
	# minp0=minp0w=0
  
	# # str(Results)
	# # clusterSetupRNGstream(cl,seed=rep(seed,6))
	# # clusterSetupRNG(cl)
	# #--------------------------#
	# # the func to do permutation.
	# gamma_statistics_return <- function (b=1) {
		# u.null<-uNull[b,]	# extract a row here.
		# if (gamma[g]<Inf) {
		  # T1s<-sum(u.null^gamma[g])
		  # T1sw<-sum((u.null%*% with(svd(ginv(a) ),u %*% diag(sqrt(d)) ))^gamma[g])
		# } else if (gamma[g]==Inf) {
		  # T1s<-max(abs(u.null) );
		  # T1sw<-max(abs(u.null%*% with(svd(ginv(a) ),u %*% diag(sqrt(d)) )))
		# }
		  
		# if (gamma[g]==1){	
		  # t0<-max(u.null^2/a); 
		  # # count <- as.integerI(t0>T)
		# } else t0 <- NA	# for uniform umin/minP.
	   # return(c(b = b, T1s = T1s, T1sw = T1sw,t0 = t0))
	# }
	# # and need the matching uNull.
	# set.seed(seed) # to ensure the same samples are drawn for each gamma. i.e. set same seed on each slave. To set a seed on a series of slaves please use clusterSetupRNG. So batch result will be the same but not within each slave.
	# uNull <- mvrnorm(B,rep(0,length(U)),V)	# simulate uNull should come from multivariate N(0,V), only right for common variants.
	
  
  # # start computation for each gamma.
	# # Results <- foreach (g = rep(1,9),.packages="MASS",.verbose=FALSE) %dopar% {	# # for debug use all same g.
	# #######################################################
	# ##############start computation for each gamma#########
	# #######################################################
	# #-----------------------------------------------------#
	# Results <- foreach (g = 1:length(gamma),.packages="MASS",.verbose=FALSE,.errorhandling="pass") %dopar% {	# now this loop level are serial for purpose.
  # # for (g in 1:length(gamma)) {	# now this loop level are serial for purpose.
		# # collect_garbage()
		# # catt(g,'started!')
		# # set.seed(seed) # to ensure the same samples are drawn for each gamma. i.e. set same seed on each slave. To set a seed on a series of slaves please use clusterSetupRNG. So batch result will be the same but not within each slave.
		# results <-simplify2array( lapply(1:B,gamma_statistics_return) )
		# #
		# T_matrix_temp <- results	# for one gamma value, contains B permutations.
		# # names(gamma_statistics_return) <- gamma
		# p= mean( abs(spu[g])< abs(T_matrix_temp["T1s",]) )
		# pw= mean( abs(spuw[g])<abs(T_matrix_temp["T1sw",]) )
		
		# #
		# if (g==1) t0_vec <-T_matrix_temp["t0",]  else t0_vec <- NULL
		# #
		# P0s=(B-rank(abs(T_matrix_temp["T1s",]))) / (B-1)
		# P0sw=(B-rank(abs(T_matrix_temp["T1sw",]))) / (B-1)	
		# # return(list(minp0 = minp0, minp0w = minp0w, p = p, pw = pw ))
		# # catt(g,'finished!')
		# return(list(t0_vec = t0_vec,P0s = P0s, P0sw = P0sw, p = p, pw = pw ))
	# }	# end of foreach  
	# names(Results) <- gamma	# this obj may not need.
  # #	#######################################################
  # #	#######################################################
  # #	#######################################################

  	# #---------------------------------------------------#
	# foreach( g = 1:length(Results)) %do% {	# avoid parallel here, too easy task will cost more time on communication then computation.
		# temp <- Results[[g]]
		# if (g==1) minp0=temp$P0s else minp0[which(minp0>temp$P0s)] <- temp$P0s[which(minp0>temp$P0s)]
		# if (g==1) minp0w=temp$P0sw else minp0w[which(minp0w>temp$P0sw)] <- temp$P0sw[which(minp0w>temp$P0sw)]
		# cat("g=",g,"finished!\n")
	# }
	# p = sapply(Results, function(x) x$p)
	# pw = sapply(Results, function(x) x$pw)
	# Paspu <- mean( minp0<min(p) )
	# Paspuw <- mean( minp0w<min(pw) )
	# #
	# p_vector_output <- c(p,aspu = Paspu)
	# pw_vector_output <- c(pw,aspu.w = Paspuw)
	# #
	# umin= mean( Results[[1]][["t0_vec"]]>T )

	# return(list(spu=p_vector_output,spu.w=pw_vector_output,minP=umin))
# }


#########################################
# #=====================================#
# the older version below.
# SPU3<-function(U,V,gamma=c(1:8,Inf),B){
	# seed=12345
	# a=ifelse(diag(V)>1e-15,diag(V),1e-15)

	# ##observed spu
	# spu=spuw=rep(NA,length(gamma))
	# for (g in 1:length(gamma)){
	# if (gamma[g]<Inf) spu[g]<-sum(U^gamma[g]) else spu[g]<-max(abs(U))
	# if (gamma[g]<Inf) spuw[g]<-sum((U/sqrt(a))^gamma[g]) else spuw[g]<-max(abs(U/sqrt(a)))
	# }
	# ##UminP
	# T<-max(U^2/a)
	# cat("statistic calculated","\n")
	# ## perm 
	# p=pw=rep(NA,length(gamma))
	# # T1s=T1sw=numeric(B)
	# # count=0
	# gamma_statistics_return <- list()
	# minp0=minp0w=0
  
	# # str(Results)
	# # clusterSetupRNGstream(cl,seed=rep(seed,6))
	# # clusterSetupRNG(cl)
	# #--------------------------#
	# # the func to do permutation.
	# gamma_statistics_return <- function (b=1) {
		# u.null<-mvrnorm(1,rep(0,length(U)),V)	# sampling here.
		# if (gamma[g]<Inf) {
		  # T1s<-sum(u.null^gamma[g])
		  # T1sw<-sum((u.null/sqrt(a))^gamma[g])
		# } else if (gamma[g]==Inf) {
		  # T1s<-max(abs(u.null) );
		  # T1sw<-max(abs(u.null/sqrt(a)))
		# }
		  
		# if (gamma[g]==1){	# extra computation for gamma == 1. equivalent to minP method.
		  # t0<-max(u.null^2/a); 
		  # # count <- as.integerI(t0>T)
		# } else t0 <- NA	# for uniform umin/minP.
	   # return(c(b = b, T1s = T1s, T1sw = T1sw,t0 = t0))
	# }
  
  # # start computation for each gamma.
	# # Results <- foreach (g = rep(1,9),.packages="MASS",.verbose=FALSE) %dopar% {	# # for debug use all same g.
	# Results <- foreach (g = 1:length(gamma),.packages="MASS",.verbose=FALSE,.errorhandling="pass") %dopar% {	# now this loop level are serial for purpose.
  # # for (g in 1:length(gamma)) {	# now this loop level are serial for purpose.
		# # collect_garbage()
		# # catt(g,'started!')
		# set.seed(seed) # to ensure the same samples are drawn for each gamma. i.e. set same seed on each slave. To set a seed on a series of slaves please use clusterSetupRNG. So batch result will be the same but not within each slave.
		# results <-simplify2array( lapply(1:B,gamma_statistics_return) )
		# #
		# T_matrix_temp <- results	# for one gamma value, contains B permutations.
		# # names(gamma_statistics_return) <- gamma
		# p= mean( abs(spu[g])< abs(T_matrix_temp["T1s",]) )
		# pw= mean( abs(spuw[g])<abs(T_matrix_temp["T1sw",]) )
		
		# #
		# if (g==1) t0_vec <-T_matrix_temp["t0",]  else t0_vec <- NULL
		# #
		# P0s=(B-rank(abs(T_matrix_temp["T1s",]))) / (B-1)
		# P0sw=(B-rank(abs(T_matrix_temp["T1sw",]))) / (B-1)	
		# # return(list(minp0 = minp0, minp0w = minp0w, p = p, pw = pw ))
		# # catt(g,'finished!')
		# return(list(t0_vec = t0_vec,P0s = P0s, P0sw = P0sw, p = p, pw = pw ))
	# }	# end of foreach  
	# names(Results) <- gamma	# this obj may not need.
  # #
  	# #
	# foreach( g = 1:length(Results)) %do% {
		# temp <- Results[[g]]
		# if (g==1) minp0=temp$P0s else minp0[which(minp0>temp$P0s)] <- temp$P0s[which(minp0>temp$P0s)]
		# if (g==1) minp0w=temp$P0sw else minp0w[which(minp0w>temp$P0sw)] <- temp$P0sw[which(minp0w>temp$P0sw)]
		# cat("g=",g,"finished!\n")
	# }
	# p = sapply(Results, function(x) x$p)
	# pw = sapply(Results, function(x) x$pw)
	# Paspu <- mean( minp0<min(p) )
	# Paspuw <- mean( minp0w<min(pw) )
	# #
	# p_vector_output <- c(p,aspu = Paspu)
	# pw_vector_output <- c(pw,aspu.w = Paspuw)
	# #
	# umin= mean( Results[[1]][["t0_vec"]]>T )

	# return(list(spu=p_vector_output,spu.w=pw_vector_output,minP=umin))
# }


###############################################################################
# more comments are within comments below

# SPU3<-function(U,V,gamma=c(1:8,Inf),B){
	# seed=12345
	# a=ifelse(diag(V)>1e-15,diag(V),1e-15)

	# ##observed spu
	# spu=spuw=rep(NA,length(gamma))
	# for (g in 1:length(gamma)){
	# if (gamma[g]<Inf) spu[g]<-sum(U^gamma[g]) else spu[g]<-max(abs(U))
	# if (gamma[g]<Inf) spuw[g]<-sum((U/sqrt(a))^gamma[g]) else spuw[g]<-max(abs(U/sqrt(a)))
	# }
	# ##UminP
	# T<-max(U^2/a)
	# cat("statistic calculated","\n")
	# ## perm 
	# p=pw=rep(NA,length(gamma))
	# # T1s=T1sw=numeric(B)
	# # count=0
	# gamma_statistics_return <- list()
	# minp0=minp0w=0
  
	# # str(Results)
	# # clusterSetupRNGstream(cl,seed=rep(seed,6))
	# # clusterSetupRNG(cl)
	# #--------------------------#
	# # the func to do permutation.
	# gamma_statistics_return <- function (b=1) {
		# u.null<-mvrnorm(1,rep(0,length(U)),V)	# sampling here.
		# if (gamma[g]<Inf) {
		  # T1s<-sum(u.null^gamma[g])
		  # T1sw<-sum((u.null/sqrt(a))^gamma[g])
		# } else if (gamma[g]==Inf) {
		  # T1s<-max(abs(u.null) );
		  # T1sw<-max(abs(u.null/sqrt(a)))
		# }
		  
		# if (gamma[g]==1){	# extra computation for gamma == 1. equivalent to minP method.
		  # t0<-max(u.null^2/a); 
		  # # count <- as.integerI(t0>T)
		# } else t0 <- NA	# for uniform umin/minP.
	   # return(c(b = b, T1s = T1s, T1sw = T1sw,t0 = t0))
	# }
  
  # Results <- foreach (g = rep(1,9),.packages="MASS",.verbose=FALSE) %dopar% {	# now this loop level are serial for purpose.
  # # for (g in 1:length(gamma)) {	# now this loop level are serial for purpose.
    # set.seed(seed) # to ensure the same samples are drawn for each gamma
	# # rnorm(10)
	# # }
	
	# # clusterSetupRNGstream(cl,seed=rep(123,6)) # a vector of 6 to guarantee each slave has the same seed settings.
	# # clusterSetupRNG(cl) # a vector of 6 to guarantee each slave has the same seed settings.
	# # set.seed(seed)	# debug when below use %do%
	# # gamma_statistics_return[[g]] <- foreach (b = 1:B,.combine = rbind,.packages="MASS",.verbose=FALSE) %dopar% {
        # # u.null<-mvrnorm(1,rep(0,length(U)),V)	# sampling here.
        # # if (gamma[g]<Inf) {
          # # T1s<-sum(u.null^gamma[g])
          # # T1sw<-sum((u.null/sqrt(a))^gamma[g])
        # # } else if (gamma[g]==Inf) {
		  # # T1s<-max(abs(u.null) );
		  # # T1sw<-max(abs(u.null/sqrt(a)))
		# # }
          
        # # if (gamma[g]==1){	# extra computation for gamma == 1. equivalent to minP method.
		  # # t0<-max(u.null^2/a); 
		  # # # count <- as.integerI(t0>T)
	    # # } else t0 <- NA	# for uniform umin/minP.
	   # # return(c(T1s = T1s, T1sw = T1sw,t0 = t0))
    # # }
	# # #
	# # T_matrix_temp <- gamma_statistics_return[[g]]
	# # # names(gamma_statistics_return) <- gamma
    # # p[g] = mean( abs(spu[g])< abs(T_matrix_temp[,"T1s"]) )
    # # pw[g] = mean( abs(spuw[g])<abs(T_matrix_temp[,"T1sw"]) )
	
	# # #
    # # P0s=(B-rank(abs(T_matrix_temp[,"T1s"]))) / (B-1)
    # # P0sw=(B-rank(abs(T_matrix_temp[,"T1sw"]))) / (B-1)
	
	# # #
    # # if (g==1) minp0=P0s else minp0[which(minp0>P0s)] <- P0s[which(minp0>P0s)]
    # # if (g==1) minp0w=P0sw else minp0w[which(minp0w>P0sw)] <- P0sw[which(minp0w>P0sw)]
    # # cat("g=",g,"finished!\n")
	# # # return()	# return are already in p[], pw[] and minp0/w
	# # #############
	# # ## original code
	# # # p[g]=sum(abs(spu[g])<abs(T1s))/B
    # # # pw[g]=sum(abs(spuw[g])<abs(T1sw))/B
    # # # P0s=(B-rank(abs(T1s)))/(B-1)
    # # # P0sw=(B-rank(abs(T1sw)))/(B-1)
    # # # if (g==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
    # # # if (g==1) minp0w=P0sw else minp0w[which(minp0w>P0sw)]=P0sw[which(minp0w>P0sw)]
    # # # cat("g=",g,"\n")
	# # ##
	# # #############
  # # }

	
	# # execute
	# results <-simplify2array( (1:B,gamma_statistics_return) )
	# # clusterEvalQ(cl,library(MASS))
	# # system.time( results1 <- times(60000) %dopar% gamma_statistics_return() )
	   # # # user  system elapsed 
 # # # 24.429   0.396  24.970 
  # # Timing stopped at: 0.036 0 241.013	# too slow
 	# # #---------------------------#
	# # system.time( results2 <- clusterApply(cl,1:60000,gamma_statistics_return) )
	   # # user  system elapsed  
 # # 62.204  12.873  86.513  
 # # Timing stopped at: 0 0 633.221	# too slow
 	# # #---------------------------#
	# # system.time( results3 <- clusterApplyLB(cl,1:60000,gamma_statistics_return) )
	   # # user  system elapsed
# # 185.364  25.285 220.175
# # Timing stopped at: 294.57 36.243 595.603	# too slow
	# # #---------------------------#
	# # set.seed(123)	# verified, even given same seed, mcXXX func produced different results.
	# # system.time(results4 <-mclapply(1:60000,gamma_statistics_return) )
	 # # user  system elapsed 
# # 4.676   8.996   4.552
  # # user  system elapsed 
# # 57.864   9.214   7.840 
	# # #---------------------------#
	# # set.seed(123) 	# verified, even given same seed, mcXXX func produced different results.
	# # system.time(results5 <- simplify2array( mcmapply(gamma_statistics_return,1:60000,SIMPLIFY=F) ) )
	  # # user  system elapsed
 # # 0.636   1.528   3.556
    # # user  system elapsed	# 6w, 2 cores
  # # 5.256   2.132  17.223
     # # user  system elapsed	# 6w, 30 cores
 # # 23.953  44.047  16.913
	# # #---------------------------#
	# # system.time( results6 <-lapply(1:60000,gamma_statistics_return))
	   # # user  system elapsed
  # # 4.588   0.000   4.586
     # # user  system elapsed   
 # # 29.366   0.012  29.374 
	# # #---------------------------#
	# # system.time( results7 <-foreach (b = 1:60000,.combine = rbind,.packages="MASS",.verbose=FALSE) %do% {gamma_statistics_return() } )
   # # user  system elapsed 
 # # 16.853   0.000  16.853 
    # # user  system elapsed
# # 108.363   0.004 108.362
	# # #---------------------------#
	# # system.time( results8 <-foreach (b = 1:60000,.combine = rbind,.packages="MASS",.verbose=FALSE) %dopar% {gamma_statistics_return() } )
	   # # user  system elapsed
 # # 27.538   0.448  28.136
    # # user  system elapsed 
# # 177.015   2.260 179.300 
	# # #---------------------------#
	# # results9 <- c()
	# # system.time( for (b in 1:60000) {results9[[b]] <- gamma_statistics_return(b) } )
	   # # user  system elapsed	# 10000
  # # 5.704   0.040   5.742
     # # user  system elapsed 	# 60000
 # # 48.984   0.100  49.080 	
	# # #---------------------------#
	# # system.time( results10 <- parLapply(cl,1:60000,gamma_statistics_return) )
	   # # user  system elapsed	# 6w
 # # 17.041   2.568  46.739
	# # #---------------------------#
	# # system.time( results11 <- parLapplyLB(cl,1:60000,gamma_statistics_return) )
   # # user  system elapsed	# 6w
 # # 20.918   2.984  79.047
	
	# # #---------------------------#
	# # system.time( results12 <- clusterMap(cl,gamma_statistics_return,1:60000,.scheduling = 'static') )
 # # Timing stopped at: 0.004 0 148.402	# too slow.
	
	# # #---------------------------#
	# # system.time( results13 <- Map(gamma_statistics_return,1:60000) )
   # # user  system elapsed 
 # # 30.282   0.004  30.281 
	
		# # #---------------------------#
	# # system.time( results14 <- pvec(1:60000,gamma_statistics_return) )
	# # # error since must be # of tasks = # of cores assigned.
		
	# #
	# T_matrix_temp <- results	# for one gamma value, contains B permutations.
	# # names(gamma_statistics_return) <- gamma
    # p= mean( abs(spu[g])< abs(T_matrix_temp["T1s",]) )
    # pw= mean( abs(spuw[g])<abs(T_matrix_temp["T1sw",]) )
	
	# #
    # P0s=(B-rank(abs(T_matrix_temp["T1s",]))) / (B-1)
    # P0sw=(B-rank(abs(T_matrix_temp["T1sw",]))) / (B-1)
	
	
	# # return()	# return are already in p[], pw[] and minp0/w
	# #############
	# ## original code
	# # p[g]=sum(abs(spu[g])<abs(T1s))/B
    # # pw[g]=sum(abs(spuw[g])<abs(T1sw))/B
    # # P0s=(B-rank(abs(T1s)))/(B-1)
    # # P0sw=(B-rank(abs(T1sw)))/(B-1)
    # # if (g==1) minp0=P0s else minp0[which(minp0>P0s)]=P0s[which(minp0>P0s)]
    # # if (g==1) minp0w=P0sw else minp0w[which(minp0w>P0sw)]=P0sw[which(minp0w>P0sw)]
    # # cat("g=",g,"\n")
	# ##
	# #############
	# # return(list(minp0 = minp0, minp0w = minp0w, p = p, pw = pw ))
	# return(list(P0s = P0s, P0sw = P0sw, p = p, pw = pw ))
  # }	# end of foreach  
	# names(Results) <- gamma	# this obj may not need.
  # #
  	# #
	# if (g==1) minp0=P0s else minp0[which(minp0>P0s)] <- P0s[which(minp0>P0s)]
	# if (g==1) minp0w=P0sw else minp0w[which(minp0w>P0sw)] <- P0sw[which(minp0w>P0sw)]
	# cat("g=",g,"finished!\n")

	# cat("P0s caculated","\n")
	# Paspu <- mean( minp0<min(p) )
	# Paspuw <- mean( minp0w<min(pw) )
	# #
	# p_vector_output <- c(p,Paspu)
	# pw_vector_output <- c(pw,Paspuw)
	# #
	# umin= mean( Results[[1]]["t0",]>T )

	# return(list(spu=p_vector_output,spu.w=pw_vector_output,minP=umin))
# }
# ###############################################################################
