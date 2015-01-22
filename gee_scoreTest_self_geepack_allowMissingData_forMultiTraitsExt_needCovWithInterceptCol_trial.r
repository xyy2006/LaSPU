# updated on 11/25/2014 to update score test within 'geescore_for' to be more
# robust and more powerful.
library(MASS)
library(geepack)
library(gdata)
library(plyr)


# ------------the score from Test5v1, more stable with extreme condition
# cycle-------------------# ------------considered the singular value in the
# CovS, so it's not full rank, the pTscore adjust for this, thus more powerful
# than naive score test with just ginv(CovS)-------------------# score test: will
# be used in 'geescore_for'.
Score_robust <- function(U, CovS) {
    if (is.null(dim(CovS))) {
        # only one-dim:
        Tscore <- sum(U^2/CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) 
            Tscore <- 0
        pTscore <- as.numeric(1 - pchisq(Tscore, 1))
        return(data.frame(df = 1, Score.statistics = Tscore, pval = pTscore))
    } else {
        ## gInv of CovS:
        CovS.edecomp <- eigen(CovS)
        CovS.rank <- sum(abs(CovS.edecomp$values) > 1e-08)
        inveigen <- ifelse(abs(CovS.edecomp$values) > 1e-08, 1/CovS.edecomp$values, 
            0)
        P <- solve(CovS.edecomp$vectors)
        gInv.CovS <- t(P) %*% diag(inveigen) %*% P
        Tscore <- t(U) %*% gInv.CovS %*% U
        pTscore <- as.numeric(1 - pchisq(Tscore, CovS.rank))
        return(data.frame(df = CovS.rank, Score.statistics = Tscore, pval = pTscore))
    }
}



# Note: this version differ to
# gee_scoreTest_self_geepack_allowMissingData_forMultiTraitsExt_needCovWithInterceptCol.r
# in its adoption of foreach loop to do UI summation instead of faster but
# memory-consuming dlpy.
geescore_for <- function(nCluster = NULL, data, id, y, cov = "intercept", par, waves, 
    corstr = c("independence", "exchangeable", "ar1", "unstructured"), family = c("binomial", 
        "gaussian"), allowCovNA = TRUE, allowParNA = TRUE, allowYNA = TRUE, biasCorrection = T, 
    commonVariance = T, .parallel = FALSE) {
    ######################### obsolete arguments: allowMisspecified=T # doesn't matter any more.
    
    # In geepack, corstr allow: a character string specifying the correlation
    # structure.  The following are permitted: 'independence', 'exchangeable', 'ar1',
    # 'unstructured', 'userdefined', and 'fixed' browser()
    if (!.parallel) 
        `%dopar%` <- `%do%`  # change within function, for safty, ensure within func, all serial.
    stopifnot(is.data.frame(data))  # we need data to be data.frame instead of matrix.
    data <- tbl_dt(data)  # transform to data.table for speed.
    # data <- na.omit(data) # the geese/geeglm program doesn't allow NA in y or cov
    # part either, so if there is NA, delete rowwisely.  if in y, the fit and
    # residual will be N - n, where n is the # of missing rows. So ee matrix need to
    # be remodified.  if in cov, geese/geeglm will report error whatever na.action
    # equates.  if in par, already solved.
    if (!is.numeric(data[, id, with = F])) {
        firstWavesSeq <- data[, list(N = .N), by = id][, N]
        setnames(data, id, sprintf("%s_originalBackup", id))  # rename invalid id.
        newID <- rep(1:length(firstWavesSeq), times = firstWavesSeq)
        # data[,eval(id):=newID] # assign new id.
        data[[eval(id)]] <- newID  # assign new id.
    }
    if (id == "id") {
        data <- plyr::arrange(data, id, eval(as.name(waves)))  # plyr allow expr in ...
    } else data <- plyr::arrange(data, eval(as.name(id)), eval(as.name(waves)))
    # ---- backup optional cmd---# arrange_call <- (sprintf('plyr::arrange(data, %s,
    # %s)', id, waves)) data <- eval(parse(text = arrange_call))
    
    # -----------#
    if (allowYNA) {
        y_NA_index <- with(data, is.na(eval(as.name(y))))  # y is always a single col, as required in long format.
        data <- data[!y_NA_index, ]  # only filger NAs in y_NA_index.
    } else tryCatch(na.fail(data[, y, with = F]), error = stop("NA in y(response) are not allowed! Please either enable 'allowYNA' or filter/impute 'y' part by yourself!"), 
        warning = identity)
    # ---# if (cov[1] != 1) { # in multipleTraits scenarios, we dont need this.
    if (allowCovNA) {
        # cov_NA_index <- with(data, is.na(eval(as.name(cov))) )
        cov_NA_index <- apply(data[, cov, with = F], 1, function(x) any(is.na(x)))  # for data.table.
        data <- data[!cov_NA_index, ]
    } else tryCatch(na.fail(data[, cov, with = F]), error = stop("NA in cov(covariates) are not allowed! Please either enable 'allowCovNA' or filter/impute 'cov' part by yourself!"), 
        warning = identity)
    # } ---# par_NA_index <- with(data, is.na(eval(as.name(par))) ) # single par,
    # which possibly not gonna happen all the time.  par_NA_index <- with(data,
    # is.na(eval(as.name(par))) )
    par_NA_index <- apply(data[, par, with = F], 1, function(x) any(is.na(x)))  # for data.table.
    if (allowParNA & sum(par_NA_index) != 0) {
        # we didn't filter any row due to Par NAs.
        message("There are NAs in the 'par' part, with 'allowParNA = TRUE', we will keep going as usual! (results are still reliable given small amount of NAs, it's your responsibility to evaluate how many NAs there are in 'par' part. Impute your 'par' part (usually they are genetic variants) if there are too many missings to keep the accuracy of the program.)")
    } else if (!allowParNA & sum(par_NA_index) != 0) 
        {
            stop("NA in par(parameters of interests, usually genetic variants) are not allowed! Please either enable 'allowParNA' or impute 'par' part by yourself!")
        }  # if sum(par_NA_index)==0, allowParNA doesn't matter any longer.
    # we just go ahead to below part.
    
    # ==============================#
    family <- match.arg(family)
    corstr <- match.arg(corstr)
    if (is.null(nCluster)) 
        nCluster <- length(unique(data[[eval(id)]]))  # for data.frame data, how many nCluster we have after all data filtering.
    
    ####################### ---------------------# this part not really used in the below code, but can be
    ####################### useful. So keep it here.
    waves_by_id <- table(data[[eval(id)]])
    if (all(waves_by_id == max(waves_by_id))) {
        full_obs_flag <- TRUE  # perfect, we have full obs data, which will facilitate the below computation.
    } else full_obs_flag <- FALSE
    
    
    # ---------------------#
    
    
    ####################### ---------------------# this is used below for building Rw.
    different_waves_we_have <- sort(unique(waves_by_id))  # e.g. 4 or [1] 1 2 3 4 depending on full obs or missing obs occurs.
    # ---------------------#
    
    
    
    # if (cov == 1) { # only intercept there.  x <- data.frame(intercept = 1,
    # data[,par],check.names=F) } else x<- data.frame(intercept = 1,
    # data[,c(cov,par)],check.names=F) # for the latter dimension usage for whole U .
    # attr(x,'id') <- data[eval(id)] # the corresponding id for each row in x
    
    # ------------fit the model under null--------------#
    form <- as.formula(paste(y, "~ 0 +", paste(cov, collapse = "+"), sep = ""))  # y is 'char'
    # cov <- c(1,'X.2') # for debug purpose.  if (binary==T) geefit_result <- geeglm(
    # form, id = eval(id),waves = eval( as.name(waves) ),data = data,
    # corstr=corstr,family='binomial') else geefit_result <- geeglm( form, id =
    # eval(id),waves = eval( as.name(waves) ),data = data,
    # corstr=corstr,family='gaussian') # which also contains the $geese for above.
    # if (binary==T) geefit_result <- geeglm( form, id = eval( as.name(id) ),waves =
    # eval( as.name(waves) ),data = data, corstr=corstr,family='binomial') else
    # geefit_result <- geeglm( form, id = eval( as.name(id) ),waves = eval(
    # as.name(waves) ),data = data, corstr=corstr,family='gaussian') # which also
    # contains the $geese for above.
    if (id == "id") {
        geefit_result <- geeglm(form, id = id, waves = eval(as.name(waves)), data = data, 
            corstr = corstr, family = family)  # which also contains the $geese for above.
    } else geefit_result <- geeglm(form, id = eval(as.name(id)), waves = eval(as.name(waves)), 
        data = data, corstr = corstr, family = family)  # which also contains the $geese for above.
    fit0 <- geefit_result  # that is the model fitting under NULL hypothesis.
    # scall <- match.call() mnames <- c('', '', 'formula', 'data', 'offset',
    # 'weights', 'subset', 'na.action', 'id', 'waves', 'corp') cnames <- names(scall)
    # cnames <- cnames[match(mnames, cnames, 0)] mcall <- scall[cnames] if
    # (is.null(mcall$id)) mcall$id <- as.name('id') mcall[[1]] <-
    # as.name('model.frame') mcall[[2]] <- form m <- eval(mcall, parent.frame()) id
    # <- model.extract(m, id) waves <- model.extract(m, 'waves') --------extract
    # various infomation from model fit---------# R=fit0$working.correlation
    # max_longiN <- max(data[eval(waves)] ) # the max measurements within one
    # subject.  max_longiN <- data[,max(eval(as.name(waves)))] # wrong, not correct
    # for univariate trait data.
    max_longiN <- max(waves_by_id)  # the max measurements within one subject.\t# syntax for data.table obj.
    full_dim_vec <- 1:max_longiN  # for later usage in building Rws and vys.
    
    y_hat <- fitted(fit0)  # only for binomial case usage, length = n x k, e.g. 36720
    # U1<-matrix(0,ncol=1,nrow=ncol(x)) I<-matrix(0,ncol=ncol(x),nrow=ncol(x))
    e = fit0$residual  # length = n x k, a.k.a the y_resid
    data <- transform(data, y_resid = e)  # not dplyr::mutate, e is evaluate in data only. transform or plyr::mutate suffice.
    
    
    ## First get the whole score vector and information matrix under null
    if (full_obs_flag) {
        ee <- matrix(e, ncol = max_longiN, byrow = T)  # residual data frame with dim = n * ni
        # here only for below vy.  only work for non-missing data case.
    } else {
        # ---a sub func here--- # construct the ee residual matrix to n*ni.
        make_ee_row <- function(data_chunk = subset(data, subset = eval(as.name(id)) == 
            1)) {
            # browser()
            measurements_seq <- data_chunk[[waves]]
            y_resid_vector <- rep(NA, max_longiN)
            y_resid_vector[measurements_seq] <- data_chunk$y_resid
            return(y_resid_vector)
        }
        if (id == "id") {
            ee <- daply(.data = data, .variables = .(id), .fun = make_ee_row, .progress = "text", 
                .parallel = FALSE, .paropts = list(.packages = c("MASS", "gdata", 
                  "plyr")), .inform = F)  # ply family consume extreme large memory when we have large number of variants, so avoid use parallel scheme will alleviate this issue.
        } else {
            ee <- daply(.data = data, .variables = .(eval(as.name(id))), .fun = make_ee_row, 
                .progress = "text", .parallel = FALSE, .paropts = list(.packages = c("MASS", 
                  "gdata", "plyr")), .inform = F)
        }
        
    }
    # ===================================# ---a sub func here--- # collect all
    # possible measurements_seq categories in the data.  e.g. 1234, 123, 134,12,24,
    # etc.  only do it when we have missing data.
    if (!full_obs_flag) {
        collect_all_measurements_seqs <- function(data_chunk = subset(data, subset = eval(as.name(id)) == 
            1)) {
            # browser()
            measurements_seq <- data_chunk[[waves]]
            return(measurements_seq)
        }
        if (id == "id") {
            measurements_seq_pool <- dlply(.data = data, .variables = .(id), .fun = collect_all_measurements_seqs, 
                .progress = "text", .parallel = FALSE, .paropts = list(.packages = c("MASS", 
                  "gdata", "plyr")), .inform = F)
        } else {
            measurements_seq_pool <- dlply(.data = data, .variables = .(eval(as.name(id))), 
                .fun = collect_all_measurements_seqs, .progress = "text", .parallel = FALSE, 
                .paropts = list(.packages = c("MASS", "gdata", "plyr")), .inform = F)
        }
        
        # measurements_seq_pool # will be used later in building subject specific Rw.
        measurements_seq_pool_pasted <- sapply(measurements_seq_pool, paste, collapse = "")  # for N subjects, length = n
        measurements_seq_categories <- names(table(measurements_seq_pool_pasted))
    }
    
    if (commonVariance & family == "gaussian") {
        # v-cov y is exactly v-cov e vy_full<-cov(ee,use='pairwise.complete.obs') # k by
        # k v-cov matrix; pairwise.complete.obs retains the most informations. k =
        # max_longiN vy_full<- matrix( rowSums(apply(ee,1,tcrossprod),na.rm=T)/
        # (nCluster),nrow = max_longiN,byrow=T) # k by k v-cov matrix; commonVariance
        # only for quantitative case, binomial cannot adopt this. # biased, see this
        # example show: > apply(aaa,1,tcrossprod) [,1] [,2] [,3] [1,] 1 4 9 [2,] 4 NA 18
        # [3,] 7 16 27 [4,] 4 NA 18 [5,] 16 NA 36 [6,] 28 NA 54 [7,] 7 16 27 [8,] 28 NA
        # 54 [9,] 49 64 81 > rowSums(apply(aaa,1,tcrossprod)) [1] 14 NA 50 NA NA NA 50 NA
        # 194 > rowMeans(apply(aaa,1,tcrossprod)) [1] 4.666667 NA 16.666667 NA NA NA
        # 16.666667 [8] NA 64.666667 > rowMeans(apply(aaa,1,tcrossprod), na.rm=T) [1]
        # 4.666667 11.000000 16.666667 11.000000 26.000000 41.000000 16.666667 [8]
        # 41.000000 64.666667 > rowSums(apply(aaa,1,tcrossprod), na.rm=T)/ 3 [1] 4.666667
        # 7.333333 16.666667 7.333333 17.333333 27.333333 16.666667 [8] 27.333333
        # 64.666667 or cross-sectional data.
        if (max_longiN == 1) {
            vy_full <- mean(ee^2, na.rm = T)  # scalar
            # Browse[2]> identical(as.vector(ee^2), apply(ee,1,tcrossprod)) [1] TRUE
            
        } else vy_full <- matrix(rowMeans(apply(ee, 1, tcrossprod), na.rm = T), nrow = max_longiN, 
            byrow = T)  # k by k v-cov matrix; commonVariance only for quantitative case, binomial cannot adopt this.
        # ---# get various vy_partial for missing obs scenarios.
        if (!full_obs_flag) {
            if (max_longiN == 1) {
                # or cross-sectional data.
                vys <- vy_full
            } else {
                vys <- foreach(msc = measurements_seq_categories) %do% {
                  # no need to parallel if (nchar(msc) > 1) { # dont need this, as even one obs in
                  # either longitudinal missing data setting or cross-sectional setting, it has a
                  # vy comming from common vys.
                  dim_vec <- as.integer(strsplit(msc, split = "")[[1]])
                  dim_select_index <- full_dim_vec %in% dim_vec
                  vy = vy_full[dim_select_index, dim_select_index]
                  # } else vy = NULL
                  return(vy)
                }
                names(vys) <- measurements_seq_categories  # now we have the Rws storing all needed Rw types for data.
            }
        }
        # this cov func is the same as below. I.e., Sum_i { (ee[i, ] - ee_mean) %*%
        # t(ee[i, ] - ee_mean) } / (N - 1). N is the nCluster and also the nrow(ee).
        # ee_mean <- colMeans(ee) ee_centered <- sweep(ee,2,ee_mean) matrix(
        # rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/ (nCluster-1 ),nrow =
        # max_longiN,byrow=T) # the estimates of cov(r), r is the matrix of ri in a row.
        # Browse[2]> matrix( rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/
        # (nCluster-1 ),nrow = max_longiN,byrow=T) [,1] [,2] [,3] [,4] [1,] 2.000849
        # 1.731253 1.476511 1.328423 [2,] 1.731253 1.978508 1.627722 1.438636 [3,]
        # 1.476511 1.627722 1.854880 1.620617 [4,] 1.328423 1.438636 1.620617 1.984064
        
        # Browse[2]> vy [,1] [,2] [,3] [,4] [1,] 2.000849 1.731253 1.476511 1.328423 [2,]
        # 1.731253 1.978508 1.627722 1.438636 [3,] 1.476511 1.627722 1.854880 1.620617
        # [4,] 1.328423 1.438636 1.620617 1.984064
    } else {
        vy_full <- NULL
        # ee_mean <- colMeans(ee,na.rm=T) # for later on individual i use to substract. #
        # dont need this.
    }
    
    # if (binary==T) v_of_mu=y_hat*(1-y_hat) else v_of_mu=rep(1,nrow(data)) # v_of_mu
    # is the A in the formula, where Vi = Ai(1/2) Rwi Ai(1/2)
    v_of_mu <- switch(family, gaussian = rep(1, nrow(data)), binomial = y_hat * (1 - 
        y_hat), stop("Please specify a correct type of data distribution among \"gaussian\" and \"binomial\"!"))
    # data <- mutate(data,X.intercept = 1,y_hat = y_hat,v_of_mu=v_of_mu) # we let
    # data now include the col for intercept.
    data <- transform(data, y_hat = y_hat, v_of_mu = v_of_mu)  # we dont need intercept, since they are already included in each row for each trait. (heterogeneous intercepts)
    # change mutate to transform for safety reason.
    
    
    ################################################ construct the Rw################## working correlation matrix########
    ################################################ -------define the Rw for max_longiN, then reduced dimension can be subsetted
    ################################################ out.--------------# Rws <- foreach ( d_w = different_waves_we_have) %do% { # no
    ################################################ need to parallel at all.  Rws <- foreach ( d_w = max_longiN) %do% { # no need
    ################################################ to parallel at all.  if (d_w > 1) {
    dim_R <- max_longiN  # the nrow and ncol for Rw of this subject i.
    if (dim_R > 1) {
        # only when in longitudinal setting
        if (corstr == "independence") {
            rho0 <- 0
        } else rho0 <- fit0$geese$alpha
        # for cs and ar1, a scalar returned.  for unstructured, a vector returned e.g.
        # Browse[2]> geefit_result$geese$alpha alpha.1:2 alpha.1:3 alpha.1:4 alpha.2:3
        # alpha.2:4 alpha.3:4 0.8511 0.7668 0.6997 0.8819 0.7828 0.8822
        
        
        # ------building the Rw for subject i---------------# work for no-missing data
        # case.
        independence_Rw = diag(1, nrow = dim_R, ncol = dim_R)
        ar1_Rw = rho0^abs(row(independence_Rw) - col(independence_Rw))  # for ar1
        cs_Rw = rho0^(abs(row(independence_Rw) - col(independence_Rw)) > 0)  # for exchangeable
        # unstructured is a little bit more work
        unstructured_Rw = independence_Rw  # need gdata package
        lowerTriangle(unstructured_Rw) <- rho0
        for (j in 2:nrow(independence_Rw)) {
            for (i in 1:(j - 1)) {
                unstructured_Rw[i, j] <- unstructured_Rw[j, i]
            }
        }
        # now we have unstructured_Rw
        
        # -------give appropriate Rw for max_longiN------------------#
        Rw_full_dim = switch(corstr, independence = independence_Rw, ar1 = ar1_Rw, 
            exchangeable = cs_Rw, unstructured = unstructured_Rw, stop("Please provide an explicit correlation structure among \"independence\",\"exchangeable\",\"ar1\" or \"unstructured\"!"))
    } else Rw_full_dim <- NULL
    # } # else Rw = NULL return(Rw) } names(Rws) <- max_longiN # each Rw within Rws
    # has a name == different_waves_we_have for subject i.  Browse[4]> Rws $`4` [,1]
    # [,2] [,3] [,4] [1,] 1.0000000 0.8813135 0.7767136 0.6845282 [2,] 0.8813135
    # 1.0000000 0.8813135 0.7767136 [3,] 0.7767136 0.8813135 1.0000000 0.8813135 [4,]
    # 0.6845282 0.7767136 0.8813135 1.0000000
    Rws <- 1  # for cross-sectional data.
    if (max_longiN > 1 & (!full_obs_flag)) {
        
        Rws <- foreach(msc = measurements_seq_categories) %do% {
            # no need to parallel if (nchar(msc) > 1) {
            dim_vec <- as.integer(strsplit(msc, split = "")[[1]])
            dim_select_index <- full_dim_vec %in% dim_vec
            Rw = Rw_full_dim[dim_select_index, dim_select_index]
            # } else Rw = NULL
            return(Rw)
        }
        names(Rws) <- measurements_seq_categories  # now we have the Rws storing all needed Rw types for data.
    }
    # 
    
    # End: construct the Rw###########
    
    
    # #---- set up the empty container for U and I under null hypothesis.
    # U<-matrix(0,ncol=1,nrow=ncol(x)) # the U full under null #------AK is -EU(1
    # derive); BK is EUU'.----------# AK<-matrix(0,ncol=ncol(x),nrow=ncol(x))
    # BK<-matrix(0,ncol=ncol(x),nrow=ncol(x)) I <-
    # matrix(0,ncol=ncol(x),nrow=ncol(x)) # # The cov for above U, same as BK should
    # be. ncol robust for both df and matrix. length in matrix will be product of the
    # nrow*ncol.
    
    
    ######################################### ------derive Uifor each nfam number----# a func to derive the U vector for each
    ######################################### subject i.  this is for missingData part of the data.
    derive_Ui <- function(data_chunk) {
        # browser()
        token_toGetRw <- paste(data_chunk[[waves]], collapse = "")
        v = data_chunk$v_of_mu
        # x_allPart = x[attr(x,'id') ==unique(data_chunk[[id]]),] # including intercept
        # col, cov if any and par.  if (cov[1] == 1) { # only intercept there.  x_allPart
        # = subset(data_chunk,select=c(par)) # including intercept col, cov if any and
        # par. # dplyr syntax # x_allPart = data_chunk[,c(par),with=F] # data.table
        # syntax.  } else x_allPart = data_chunk[,c(cov,par),with=F]
        x_allPart = subset(data_chunk, select = c(cov, par))  # in multipleTraits, we all need this, cov at least contains intercepts.
        
        
        x_allPart <- as.matrix(x_allPart)  # a matrix, e.g. 1-4 row.
        dev <- x_allPart * v  # the D in the U function
        # x_allPart * c(.5,.4,.6,.7) == diag(c(.5,.4,.6,.7)) %*% x_allPart
        y_resid <- data_chunk$y_resid  # residual vector, i.e. yi - yhati. length = k
        # if (exists('vys') ) { # we have vys in parent.frame iff commonVariance & family
        # == 'gaussian' as specified in line.177 # too slow for exists.  we have vys in
        # parent.frame iff commonVariance & family == 'gaussian' as specified in line.177
        if (commonVariance & family == "gaussian") {
            vy = vys[[token_toGetRw]]
        } else {
            # ee_mean_sub <- ee_mean[data_chunk[[waves ]]] # the corresponding ee_mean_sub
            # part according to waves. #
            vy = y_resid %*% t(y_resid)  # if not commonVariance
        }
        
        # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame if
        # nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.
        Rw = Rws[[token_toGetRw]]  # Rws comming from parent.frame
        # use token_toGetRw to get the Rws corresponding element.
        
        # ---------------start calculate U and Sigma of U------------------------#
        # equivalent to nchar(token_toGetRw) > 1
        if (nrow(data_chunk) > 1) {
            V <- diag(sqrt(v)) %*% Rw %*% diag(sqrt(v))  # the V_i
            U0 <- t(dev) %*% ginv(V) %*% y_resid  # for one subject, U full function under the null.
            # ------AK is -EU(1 derive); BK is EUU'.----------# AK and BK not needed any
            # more, so comment to save time.  if Rw() is correctly specified, then AK = BK
            # AKi = t(dev) %*% ginv(V) %*% dev
            I0 <- t(dev) %*% ginv(V) %*% vy %*% t(ginv(V)) %*% dev  # assume common variance vy across samples. The Variance of above U0.
            # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy) %*%t(ginv(V))%*%dev # not assume common
            # variance across samples BKi = I0
            
        } else if (nrow(data_chunk) == 1) {
            # all goes in to scalar.
            V <- v^2  # scalar, since one trait, no R is needed.
            U0 <- matrix(dev, ncol = 1) * c(y_resid)/V
            I0 <- matrix(dev, ncol = 1) %*% matrix(dev, nrow = 1) * c(vy)/V/V
            # BKi <- I0 AKi <- matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V
        }
        if (!commonVariance) 
            vy <- NULL  # setting back to NULL for safty.
        # return(list(U = U0, I = I0, AK = AKi, BK = BKi) )
        return(cbind(U = U0, I = I0))
        
        
    }
    ######################################### ------derive Uifor each nfam number----# a func to derive the U vector for each
    ######################################### subject i. updated for using dt and dplyr this is for complete part of the
    ######################################### data.  derive_Ui_forFullObs <- function(data_chunk = subset(data,subset =
    ######################################### eval(as.name(id)) == 1) ){ browser()
    derive_Ui_forFullObs <- function(data_chunk) {
        # browser() token_toGetRw <- paste(data_chunk[[waves ]], collapse='')
        v = data_chunk$v_of_mu
        # v = data_chunk[,v_of_mu] # dt syntax if (cov[1] == 1) { # only intercept there.
        # x_allPart = subset(data_chunk,select=c(par)) # including intercept col, cov if
        # any and par. # dplyr syntax # x_allPart = data_chunk[,c(par),with=F] #
        # data.table syntax.  } else { x_allPart = data_chunk[,c(cov,par),with=F]
        x_allPart = subset(data_chunk, select = c(cov, par))  # at least cov has intercept.1-4
        # }
        x_allPart <- as.matrix(x_allPart)  # a matrix, e.g. 1-4 row.
        dev <- x_allPart * v  # the D in the U function
        # x_allPart * c(.5,.4,.6,.7) == diag(c(.5,.4,.6,.7)) %*% x_allPart
        y_resid <- data_chunk$y_resid  # residual vector, i.e. yi - yhati. length = k
        # y_resid <- data_chunk[,y_resid] # residual vector, i.e. yi - yhati. length = k
        # # dt syntax if (exists('vys') ) { # we have vys in parent.frame iff
        # commonVariance & family == 'gaussian' as specified in line.177 # too slow for
        # exists.  we have vys in parent.frame iff commonVariance & family == 'gaussian'
        # as specified in line.177
        if (commonVariance & family == "gaussian") {
            vy = vy_full
        } else {
            # ee_mean_sub <- ee_mean[data_chunk[[waves ]]] # the corresponding ee_mean_sub
            # part according to waves. #
            vy = y_resid %*% t(y_resid)  # if not commonVariance
        }
        
        # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame if
        # nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.
        Rw = Rw_full_dim  # Rws comming from parent.frame
        # use token_toGetRw to get the Rws corresponding element. ---------------start
        # calculate U and Sigma of U------------------------# equivalent to
        # nchar(token_toGetRw) > 1
        if (nrow(data_chunk) > 1) {
            V <- diag(sqrt(v)) %*% Rw %*% diag(sqrt(v))  # the V_i
            U0 <- t(dev) %*% ginv(V) %*% y_resid  # for one subject, U full function under the null.
            # ------AK is -EU(1 derive); BK is EUU'.----------# if Rw() is correctly
            # specified, then AK = BK AKi = t(dev) %*% ginv(V) %*% dev
            I0 <- t(dev) %*% ginv(V) %*% vy %*% t(ginv(V)) %*% dev  # assume common variance vy across samples. The Variance of above U0.
            # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy) %*%t(ginv(V))%*%dev # not assume common
            # variance across samples BKi = I0
            
        } else if (nrow(data_chunk) == 1) {
            # all goes in to scalar.
            V <- v^2  # scalar, since one trait, no R is needed.
            U0 <- matrix(dev, ncol = 1) * c(y_resid)/V
            I0 <- matrix(dev, ncol = 1) %*% matrix(dev, nrow = 1) * c(vy)/V/V
            # BKi <- I0 AKi <- matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V
        }
        if (!commonVariance) 
            vy <- NULL  # setting back to NULL for safty.
        # return(list(U = U0, I = I0, AK = AKi, BK = BKi) )
        return(cbind(U = U0, I = I0))
    }
    # ###########################################################
    
    
    # ########################################################################
    # ###########the below is original loop structure#########################
    # ###########confirmed the above func lead to same result as this#########
    # ########################################################################
    # #------derive Uifor each nfam number----# the for loop equivalent to above
    # derive_Ui data_by_id <- grouped_dt(data,quote(id)) # produce grouped_dt, for
    # do() usage.  browser() this below for loop will be faster than dlply operation
    # on subsetting data.frame.  U_results_all_subjects_collected <-
    # derive_Ui_forFullObs(index = 1:nCluster) # apply on full index if (id == 'id')
    # { U_results_all_subjects_collected <- dlply(.data = data, .variables =.(id),
    # .fun =derive_Ui_forFullObs, .progress = 'text', .parallel = .parallel, .paropts
    # = list(.packages=c('MASS','gdata','plyr')), .inform=F) # #
    # U_results_all_subjects_collected <- do(data_by_id,derive_Ui_forFullObs) # do is
    # slow, not fast at all in dplyr.  } else U_results_all_subjects_collected <-
    # dlply(.data = data, .variables =.(eval(as.name(id))), .fun
    # =derive_Ui_forFullObs, .progress = 'text', .parallel = .parallel, .paropts =
    # list(.packages=c('MASS','gdata','plyr')),.inform=F)
    if (full_obs_flag) {
        muy <- data$y_hat
        v_of_mu <- data$v_of_mu
        y_resid <- data$y_resid
        # lag=max_longiN
        x <- data[, c(cov, par), with = F]  # the matrix part only including cov and par.
        U <- matrix(0, ncol = 1, nrow = ncol(x))
        I <- matrix(0, ncol = ncol(x), nrow = ncol(x))
        # system.time( for (i in 1: (nrow(data)/max_longiN) ) {
        UI_matrixResults <- foreach(i = 1:(nrow(data)/max_longiN), .combine = naPlus, 
            .packages = c("MASS")) %dopar% {
            muy1 <- muy[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            v = v_of_mu[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            mux1 <- as.matrix(x[(max_longiN * (i - 1) + 1):(max_longiN * i), ])
            dev <- mux1 * v
            yy <- y_resid[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            
            if (max_longiN > 1) {
                V <- diag(sqrt(c(v))) %*% Rw_full_dim %*% diag(sqrt(c(v)))
                U0 <- t(dev) %*% ginv(V) %*% yy
                I0 <- t(dev) %*% ginv(V) %*% vy_full %*% t(ginv(V)) %*% dev
            }
            
            if (max_longiN == 1) {
                V <- v^2
                U0 <- matrix(dev, ncol = 1) * c(yy)/V
                I0 <- matrix(dev, ncol = 1) %*% matrix(dev, nrow = 1) * c(vy_full)/V/V
            }
            return(cbind(U0, I0))
        }
        U <- UI_matrixResults[, 1]
        I <- UI_matrixResults[, -1]
        
        # U<-naPlus(U,U0) I<-naPlus(I,I0) } # now we have U and I from data
        message("All data part finished computation!")
        
        # ######################################################
        # ######################################################
        # ######################################################
        # ######################################################
    } else {
        # ----1st, we need to build the vector with NA in for nCluster.-----# #
        # full_vector_nClusterLong <- rep(1:4,times=nCluster) insertNA <- fuction(x =
        # waves_by_id[[1]] ){ data[[waves]][] full_dim_vec } as we dont need to put the
        # whole data into dlply computations.
        
        id_have_missingObs_index <- as.integer(names(waves_by_id)[which(waves_by_id != 
            max_longiN)])  # the index for id have missing measurements
        # data_full_subset <- subset(data,subset = !eval(id) %in%
        # id_have_missingObs_index)
        data_full_subset <- data[!eval(as.name(id)) %in% id_have_missingObs_index, 
            ]
        # data_missingObs_subset <- subset(data,subset = eval(id) %in%
        # id_have_missingObs_index)
        data_missingObs_subset <- data[eval(as.name(id)) %in% id_have_missingObs_index, 
            ]
        # browser() d_f_subset <- data_full_subset[1:(1500*4)]
        muy <- data_full_subset$y_hat
        v_of_mu <- data_full_subset$v_of_mu
        y_resid <- data_full_subset$y_resid
        # lag=max_longiN
        x <- data_full_subset[, c(cov, par), with = F]
        U <- matrix(0, ncol = 1, nrow = ncol(x))
        I <- matrix(0, ncol = ncol(x), nrow = ncol(x))
        # system.time( for (i in 1: (nrow(data_full_subset)/max_longiN) ) {
        UI_matrixResults <- foreach(i = 1:(nrow(data_full_subset)/max_longiN), .combine = naPlus, 
            .packages = c("MASS")) %dopar% {
            muy1 <- muy[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            v = v_of_mu[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            mux1 <- as.matrix(x[(max_longiN * (i - 1) + 1):(max_longiN * i), ])
            dev <- mux1 * v
            yy <- y_resid[(max_longiN * (i - 1) + 1):(max_longiN * i)]
            
            if (max_longiN > 1) {
                V <- diag(sqrt(c(v))) %*% Rw_full_dim %*% diag(sqrt(c(v)))
                U0 <- t(dev) %*% ginv(V) %*% yy
                I0 <- t(dev) %*% ginv(V) %*% vy_full %*% t(ginv(V)) %*% dev
            }
            
            if (max_longiN == 1) {
                V <- v^2
                U0 <- matrix(dev, ncol = 1) * c(yy)/V
                I0 <- matrix(dev, ncol = 1) %*% matrix(dev, nrow = 1) * c(vy_full)/V/V
            }
            return(cbind(U0, I0))
        }
        U <- UI_matrixResults[, 1]
        I <- UI_matrixResults[, -1]
        # U<-naPlus(U,U0) # robust to NA in U I<-naPlus(I,I0) } # now we have U and I
        # from data_full_subset ) #
        message("Full data part finished computation!")
        
        # need to put before below v_of_mu, x, e reassign value.
        
        # ------we need to re-build v_of_mu, x, e (which now exclude id chunks which has
        # missing obs)--------------------# v_of_mu <-
        # data_missingObs_subset[['v_of_mu']] e <- data_missingObs_subset[['y_resid']] x
        # <- data.frame(intercept = 1, subset(data_missingObs_subset, select =
        # grep('X\\.',colnames(data_missingObs_subset),perl=T,value=T)), check.names=F)
        # # in the same format as original x.  attr(x,'id') <-
        # data_missingObs_subset[eval(id)] # the corresponding id for each row in x
        # browser() suppressWarnings( Progress disabled when using parallel plyr user
        # system elapsed 1.657 0.084 3.173 Warning messages: 1: <anonymous>: ... may be
        # used in an incorrect context: .fun(piece, ...)
        
        # 2: <anonymous>: ... may be used in an incorrect context: .fun(piece, ...)
        # ---confirmed, the above info is not any wrong in the code and process.----# so
        # use suppressWarnings to get rid of this info.  if (id == 'id') {
        # U_results_missingObs_subjects_collected_subset <- dlply(.data =
        # data_missingObs_subset, .variables =.(id), .fun =derive_Ui, .progress = 'text',
        # .parallel = .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')),
        # .inform=F) # U_results_all_subjects_collected <- dlply(.data = data, .variables
        # =.(id), .fun =derive_Ui, .progress = 'text', .parallel = .parallel, .paropts =
        # list(.packages=c('MASS','gdata','plyr')), .inform=F) } else {
        # U_results_missingObs_subjects_collected_subset <- dlply(.data =
        # data_missingObs_subset, .variables =.(eval(as.name(id))), .fun =derive_Ui,
        # .progress = 'text', .parallel = .parallel, .paropts =
        # list(.packages=c('MASS','gdata','plyr')),.inform=F) #
        # U_results_all_subjects_collected <- dlply(.data = data, .variables
        # =.(eval(as.name(id))), .fun =derive_Ui, .progress = 'text', .parallel =
        # .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')),.inform=F) } # )
        # U_results_all_subjects_collected <-
        # c(U_results_full_subjects_collected_subset,U_results_missingObs_subjects_collected_subset)
        ids <- unique(data_missingObs_subset[[id]])
        # system.time(
        message("Data with missing observation part start computation!")
        U_results_missingObs_subjects_collected_subset <- foreach(id_s = ids, .combine = naPlus, 
            .verbose = F, .errorhandling = "stop", .packages = c("MASS","dplyr")) %dopar% 
            {
                # U_results_missingObs_subjects_collected_subset <- foreach(id_s =
                # ids,.verbose=F,.errorhandling='pass') %do%{ derive_Ui(data_chunk =
                # filter(data_missingObs_subset,eval(as.name(id)) %in% id_s) )
                derive_Ui(data_chunk = filter(data_missingObs_subset, data_missingObs_subset[[id]] %in% 
                  id_s))
                # derive_Ui(data_chunk =
                # filter(data_missingObs_subset,data_missingObs_subset[[id]] %in% 32) ) %>% class
            }
        # ids <- unique(data_full_subset[[id]]) system.time(
        # U_data_full_subset_subjects_collected_subset <- foreach(id_s =
        # ids,.combine=naPlus,.verbose=F,.errorhandling='pass') %do%{ #
        # derive_Ui(data_chunk = filter(data_missingObs_subset,eval(as.name(id)) %in%
        # id_s) ) # derive_Ui_forFullObs(data_chunk =
        # filter(data_full_subset,data_full_subset[[id]] %in% id_s) )
        # derive_Ui(data_chunk = filter(data_full_subset,data_full_subset[[id]] %in%
        # id_s) ) # derive_Ui(data_chunk =
        # filter(data_missingObs_subset,data_missingObs_subset[[id]] %in% 32) ) %>% class
        # } ) # user system elapsed # for derive_Ui_forFullObs do # 364.015 0.000 363.980
        # # user system elapsed # for derive_Ui do # 364.879 0.000 364.849 #
        
        # U<-matrix(0,ncol=1,nrow=ncol(x)) I<-matrix(0,ncol=ncol(x),nrow=ncol(x)) #
        # system.time( system.time( for (i in 1: (nrow(data_full_subset)/max_longiN) ) {
        # muy1<-muy[(max_longiN*(i-1)+1):(max_longiN*i)]
        # v=v_of_mu[(max_longiN*(i-1)+1):(max_longiN*i)]
        # mux1<-as.matrix(x[(max_longiN*(i-1)+1):(max_longiN*i),]) dev<-mux1*v
        # yy<-y_resid[(max_longiN*(i-1)+1):(max_longiN*i)]
        
        # if (max_longiN>1){ V<-diag(sqrt(c(v)))%*%Rw_full_dim%*%diag(sqrt(c(v)))
        # U0<-t(dev)%*%ginv(V)%*%yy I0<-t(dev)%*%ginv(V)%*%vy_full%*%t(ginv(V))%*%dev }
        
        # if (max_longiN==1){ V<-v^2 U0<-matrix(dev,ncol=1)*c(yy)/V
        # I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy_full)/V/V }
        
        # U<-naPlus(U,U0) # robust to NA in U I<-naPlus(I,I0) } # now we have U and I
        # from data_full_subset
        
        
        # ) # user system elapsed # 205.809 0.000 209.050
        
        
        # )
        U <- naPlus(U, U_results_missingObs_subjects_collected_subset[, 1])
        I <- naPlus(I, U_results_missingObs_subjects_collected_subset[, -1])
    }
    
    # Browse[1]> names(U_results_all_subjects_collected[[1]]) [1] 'U' 'I' 'AK' 'BK'
    
    ################################################## ------very slow-------------- since large matrix manipulation in from ever I
    ################################################## step calculation U <-
    ################################################## colSums(laply(U_results_all_subjects_collected,function(x) x$U ), na.rm=T) #
    ################################################## na.rm =T is robust for missing values in par. such individual cell of missing
    ################################################## value will be skipped from summation. (cellwise skip, so if only 1 obs left,
    ################################################## the colSums will be this 1 obs value.)  I <-
    ################################################## colSums(laply(U_results_all_subjects_collected,function(x) x$I ), na.rm=T) #
    ################################################## verified the dimension is correct to sum over.  AK <-
    ################################################## colSums(laply(U_results_all_subjects_collected,function(x) x$AK ), na.rm=T)#
    ################################################## verified the dimension is correct to sum over.  BK <-
    ################################################## colSums(laply(U_results_all_subjects_collected,function(x) x$BK ), na.rm=T)#
    ################################################## verified the dimension is correct to sum over.
    
    # #---define a very useful function na+ here-------------------# naPlus <-
    # function(x,y) { x[is.na(x)] <- 0; y[is.na(x)] <- 0; x + y # this allow pairwise
    # +, not sum to a scalar.  }
    
    # #----use it on getting the U, I, AK and BK------------# U <-
    # lapply(U_results_all_subjects_collected,function(x) x$U) U <- Reduce(naPlus,U)
    
    # I <- lapply(U_results_all_subjects_collected,function(x) x$I) I <-
    # Reduce(naPlus,I)
    
    
    # ---------------------------------# since the result is the same whatever the
    # allowMisspecified is, we will not use AK and BK any more.  AK <-
    # lapply(U_results_all_subjects_collected,function(x) x$AK) AK <-
    # Reduce(`na+`,AK)
    
    # BK <- lapply(U_results_all_subjects_collected,function(x) x$BK) BK <-
    # Reduce(`na+`,BK)
    
    
    # browser() #
    
    # -------------# alpha or beta, or gammma or beta. I.e., one for X (cov) one for
    # Z (parameter).  calculate covariance matrix of gamma conditional on beta I11 is
    # the covariance matrix of beta; I22 is the covariance matrix of gamma;I12 is the
    # covariance matrix of beta-gamma
    
    # giving the right names in case in scalar case, the name will not be given
    # automatically.  if (cov[1] == 1) { covName <- 'X.intercept' } else {
    covName <- cov  # we only have one case, since intercept is already integrated in cov in multipleTraits scenarios
    # }
    U = as.matrix(U)  # to transfrom possible scalar to matrix.
    rownames(U) = c(covName, par)
    I = as.matrix(I)
    dimnames(I) = list(c(covName, par), c(covName, par))
    # AK = as.matrix(AK) BK = as.matrix(BK) dimnames(AK) = dimnames(BK) <-
    # list(c(covName,par),c(covName,par) )
    
    # ----------#
    Upar <- U[par, , drop = F]  # U of paramters (of interests) part.
    I11 <- as.matrix(I[covName, covName])  # it's also BK
    I12 <- t(as.matrix(I[par, covName]))
    I21 <- as.matrix(I[par, covName])
    I22 <- as.matrix(I[par, par])
    # if ( !allowMisspecified) { for single SNP considerations.
    if (length(par) == 1) {
        V <- I22 - I12 %*% ginv(I11) %*% I21
    } else V <- I22 - I21 %*% ginv(I11) %*% I12
    
    # } else if (allowMisspecified) { AK11 = as.matrix(AK[covName,covName]) AK12 =
    # t(as.matrix(AK[par,covName])) AK21 = as.matrix(AK[par,covName]) AK22 =
    # as.matrix(AK[par,par]) # V <- I22 - AK21 %*% ginv(AK11) %*% I12 - I21 %*%
    # ginv(AK11) %*% AK12 + AK21 %*% ginv(AK11) %*% I11 %*% ginv(AK11) %*% AK12 }
    
    # ------------execute the standard gee score test--------------# #
    # score_statistic <- t(Upar)%*%ginv(V)%*%Upar # it doesn't consider the singular
    # value in V, thus is not most powerful.  # # 38.19595 for correctly specified #
    # # 38.19595 for allow misspecified # df=length(par) #
    # pval<-1-pchisq(score_statistic,df) # the gee score test pvalue.
    # return(list(U=Upar,Cov=V,geeScoreTest_out=data.frame(df=df,
    # Score.statistics=score_statistic,pval=pval)))
    
    
    
    score_result <- Score_robust(U = Upar, CovS = V)  # it returns the df with 3 elements
    return(list(U = Upar, Cov = V, geeScoreTest_out = score_result))
}



# ------------old one without using data.table and dplyr------------------#
# geescore<-function(nCluster = NULL,data,id,y,cov =
# 1,par,waves,corstr=c('independence','exchangeable','ar1','unstructured'),family
# = c('binomial','gaussian'), allowCovNA = TRUE , allowParNA = TRUE, allowYNA =
# TRUE, allowMisspecified=T, biasCorrection=T, commonVariance=T,.parallel=FALSE){
# # In geepack, corstr allow: # a character string specifying the correlation
# structure.  The # following are permitted: 'independence', 'exchangeable', #
# 'ar1', 'unstructured', 'userdefined', and 'fixed' # browser()
# stopifnot(is.data.frame(data)) # we need data to be data.frame instead of
# matrix.  # data <- na.omit(data) # the geese/geeglm program doesn't allow NA in
# y or cov part either, so if there is NA, delete rowwisely.  # if in y, the fit
# and residual will be N - n, where n is the # of missing rows. So ee matrix need
# to be remodified.  # if in cov, geese/geeglm will report error whatever
# na.action equates.  # if in par, already solved.  if (id == 'id') { data <-
# arrange(data, id, eval(as.name(waves)) ) } else data <- arrange(data,
# eval(as.name(id)), eval(as.name(waves)) )

# #-----------# if (allowYNA){ y_NA_index <- with(data, is.na(eval(as.name(y))) )
# data <- data[!y_NA_index,] # only filger NAs in y_NA_index.  } else
# tryCatch(na.fail(data[y]), error = stop('NA in y(response) are not allowed!
# Please either enable 'allowYNA' or filter/impute 'y' part by
# yourself!'),warning = identity) #---# if (cov != 1) { # only matter when cov
# are not only intercept if (allowCovNA) { cov_NA_index <- with(data,
# is.na(eval(as.name(cov))) ) data <- data[!cov_NA_index,] } else
# tryCatch(na.fail(data[cov]), error = stop('NA in cov(covariates) are not
# allowed! Please either enable 'allowCovNA' or filter/impute 'cov' part by
# yourself!'),warning = identity) } #---# par_NA_index <- with(data,
# is.na(eval(as.name(par))) ) if (allowParNA & sum(par_NA_index)!=0 ) { # we
# didn't filter any row due to Par NAs.  message('There are NAs in the 'par'
# part, with 'allowParNA = TRUE', we will keep going as usual! (results are still
# reliable given small amount of NAs, it's your responsibility to evaluate how
# many NAs there are in 'par' part. Impute your 'par' part (usually they are
# genetic variants) if there are too many missings to keep the accuracy of the
# program.)') } else if (!allowParNA & sum(par_NA_index)!=0) { stop('NA in
# par(parameters of interests, usually genetic variants) are not allowed! Please
# either enable 'allowParNA' or impute 'par' part by yourself!') } # if
# sum(par_NA_index)==0, allowParNA doesn't matter any longer.  # we just go ahead
# to below part.

# #==============================# family <- match.arg(family) corstr <-
# match.arg(corstr) if (is.null(nCluster)) nCluster <-
# length(unique(data[[eval(id)]])) # for data.frame data, how many nCluster we
# have after all data filtering.

# ####################### #---------------------# # this part not really used in
# the below code, but can be useful. So keep it here.  waves_by_id <-
# table(data[[eval(id)]]) if (all (waves_by_id == max(waves_by_id) )) {
# full_obs_flag <- TRUE # perfect, we have full obs data, which will facilitate
# the below computation.  } else full_obs_flag <- FALSE


# #---------------------# #######################


# ####################### #---------------------# # this is used below for
# building Rw.  different_waves_we_have <- sort(unique(waves_by_id)) # e.g. 4 or
# [1] 1 2 3 4 depending on full obs or missing obs occurs.
# #---------------------# #######################



# # if (cov == 1) { # # only intercept there.  # x <- data.frame(intercept = 1,
# data[,par],check.names=F) # } else x<- data.frame(intercept = 1,
# data[,c(cov,par)],check.names=F) # for the latter dimension usage for whole U .
# # attr(x,'id') <- data[eval(id)] # the corresponding id for each row in x

# #------------fit the model under null--------------#
# form<-as.formula(paste(y,'~',paste(cov,collapse='+'),sep='') ) # y is 'char' #
# cov <- c(1,'X.2') # for debug purpose.  # if (binary==T) geefit_result <-
# geeglm( form, id = eval(id),waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family='binomial') else geefit_result <- geeglm( form, id =
# eval(id),waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family='gaussian') # which also contains the $geese for above.  #
# if (binary==T) geefit_result <- geeglm( form, id = eval( as.name(id) ),waves =
# eval( as.name(waves) ),data = data, corstr=corstr,family='binomial') else
# geefit_result <- geeglm( form, id = eval( as.name(id) ),waves = eval(
# as.name(waves) ),data = data, corstr=corstr,family='gaussian') # which also
# contains the $geese for above.  if (id == 'id') { geefit_result <- geeglm(
# form, id = id,waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family=family) # which also contains the $geese for above.  }
# else geefit_result <- geeglm( form, id = eval(as.name(id)),waves = eval(
# as.name(waves) ),data = data, corstr=corstr,family=family) # which also
# contains the $geese for above.  fit0 <- geefit_result# that is the model
# fitting under NULL hypothesis.  # scall <- match.call() # mnames <- c('', '',
# 'formula', 'data', 'offset', 'weights', 'subset', # 'na.action', 'id', 'waves',
# 'corp') # cnames <- names(scall) # cnames <- cnames[match(mnames, cnames, 0)] #
# mcall <- scall[cnames] # if (is.null(mcall$id)) # mcall$id <- as.name('id') #
# mcall[[1]] <- as.name('model.frame') # mcall[[2]] <- form # m <- eval(mcall,
# parent.frame()) # id <- model.extract(m, id) # waves <- model.extract(m,
# 'waves') #--------extract various infomation from model fit---------# #
# R=fit0$working.correlation max_longiN <- max(data[eval(waves)] ) # the max
# measurements within one subject.  full_dim_vec <- 1:max_longiN # for later
# usage in building Rws and vys.

# y_hat<-fitted(fit0) # only for binomial case usage, length = n x k, e.g. 36720
# # U1<-matrix(0,ncol=1,nrow=ncol(x)) # I<-matrix(0,ncol=ncol(x),nrow=ncol(x))
# e=fit0$residual # length = n x k, a.k.a the y_resid data <-
# data.frame(data,y_resid = e,check.names=F)


# ## First get the whole score vector and information matrix under null if
# (full_obs_flag) { ee<-matrix(e,ncol=max_longiN,byrow=T) # residual data frame
# with dim = n * ni # here only for below vy.  # only work for non-missing data
# case.  } else { # ---a sub func here--- # # construct the ee residual matrix to
# n*ni.  make_ee_row <- function(data_chunk = subset(data,subset =
# eval(as.name(id)) == 1) ){ # browser() measurements_seq <- data_chunk[[waves]]
# y_resid_vector <- rep(NA, max_longiN) y_resid_vector[measurements_seq] <-
# data_chunk$y_resid return(y_resid_vector) } if (id == 'id') { ee <- daply(.data
# = data, .variables =.(id), .fun =make_ee_row, .progress = 'text', .parallel =
# .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')), .inform=F) }
# else { ee <- daply(.data = data, .variables =.(eval(as.name(id))), .fun
# =make_ee_row, .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F) }

# } #===================================# # ---a sub func here--- # # collect all
# possible measurements_seq categories in the data.  # e.g. 1234, 123, 134,12,24,
# etc.  if (!full_obs_flag){ # only do it when we have missing data.
# collect_all_measurements_seqs <- function(data_chunk = subset(data,subset =
# eval(as.name(id)) == 1) ){ # browser() measurements_seq <- data_chunk[[waves]]
# return(measurements_seq) } if (id == 'id') { measurements_seq_pool <-
# dlply(.data = data, .variables =.(id), .fun =collect_all_measurements_seqs,
# .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')), .inform=F) } else {
# measurements_seq_pool <- dlply(.data = data, .variables =.(eval(as.name(id))),
# .fun =collect_all_measurements_seqs, .progress = 'text', .parallel = .parallel,
# .paropts = list(.packages=c('MASS','gdata','plyr')),.inform=F) }

# # measurements_seq_pool # will be used later in building subject specific Rw.
# measurements_seq_pool_pasted <- sapply(measurements_seq_pool,paste,collapse='')
# # for N subjects, length = n measurements_seq_categories <- names(
# table(measurements_seq_pool_pasted) ) }

# if (commonVariance & family == 'gaussian') { # # v-cov y is exactly v-cov e #
# vy_full<-cov(ee,use='pairwise.complete.obs') # k by k v-cov matrix;
# pairwise.complete.obs retains the most informations. k = max_longiN # vy_full<-
# matrix( rowSums(apply(ee,1,tcrossprod),na.rm=T)/ (nCluster),nrow =
# max_longiN,byrow=T) # k by k v-cov matrix; commonVariance only for quantitative
# case, binomial cannot adopt this. # biased, see this example show: # >
# apply(aaa,1,tcrossprod) # [,1] [,2] [,3] # [1,] 1 4 9 # [2,] 4 NA 18 # [3,] 7
# 16 27 # [4,] 4 NA 18 # [5,] 16 NA 36 # [6,] 28 NA 54 # [7,] 7 16 27 # [8,] 28
# NA 54 # [9,] 49 64 81 # > rowSums(apply(aaa,1,tcrossprod)) # [1] 14 NA 50 NA NA
# NA 50 NA 194 # > rowMeans(apply(aaa,1,tcrossprod)) # [1] 4.666667 NA 16.666667
# NA NA NA 16.666667 # [8] NA 64.666667 # > rowMeans(apply(aaa,1,tcrossprod),
# na.rm=T) # [1] 4.666667 11.000000 16.666667 11.000000 26.000000 41.000000
# 16.666667 # [8] 41.000000 64.666667 # > rowSums(apply(aaa,1,tcrossprod),
# na.rm=T)/ 3 # [1] 4.666667 7.333333 16.666667 7.333333 17.333333 27.333333
# 16.666667 # [8] 27.333333 64.666667 if (max_longiN == 1){ # or cross-sectional
# data.  vy_full <- mean(ee^2,na.rm=T) # scalar # Browse[2]>
# identical(as.vector(ee^2), apply(ee,1,tcrossprod)) # [1] TRUE

# } else vy_full<- matrix( rowMeans(apply(ee,1,tcrossprod),na.rm=T),nrow =
# max_longiN,byrow=T) # k by k v-cov matrix; commonVariance only for quantitative
# case, binomial cannot adopt this.  #---# # get various vy_partial for missing
# obs scenarios.  if (!full_obs_flag){ if (max_longiN == 1){ # or cross-sectional
# data.  vys <- vy_full } else { vys <- foreach ( msc =
# measurements_seq_categories ) %do% { # no need to parallel # if (nchar(msc) >
# 1) { # dont need this, as even one obs in either longitudinal missing data
# setting or cross-sectional setting, it has a vy comming from common vys.
# dim_vec <- as.integer( strsplit(msc,split='')[[1]] ) dim_select_index <-
# full_dim_vec %in% dim_vec vy = vy_full[dim_select_index,dim_select_index] # }
# else vy = NULL return(vy) } names(vys) <- measurements_seq_categories # now we
# have the Rws storing all needed Rw types for data.  } } # this cov func is the
# same as below. I.e., Sum_i { (ee[i, ] - ee_mean) %*% t(ee[i, ] - ee_mean) } /
# (N - 1). N is the nCluster and also the nrow(ee).  # ee_mean <- colMeans(ee) #
# ee_centered <- sweep(ee,2,ee_mean) # matrix(
# rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/ (nCluster-1 ),nrow =
# max_longiN,byrow=T) # the estimates of cov(r), r is the matrix of ri in a row.
# # Browse[2]> matrix( rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/
# (nCluster-1 ),nrow = max_longiN,byrow=T) # [,1] [,2] [,3] [,4] # [1,] 2.000849
# 1.731253 1.476511 1.328423 # [2,] 1.731253 1.978508 1.627722 1.438636 # [3,]
# 1.476511 1.627722 1.854880 1.620617 # [4,] 1.328423 1.438636 1.620617 1.984064

# # Browse[2]> vy # [,1] [,2] [,3] [,4] # [1,] 2.000849 1.731253 1.476511
# 1.328423 # [2,] 1.731253 1.978508 1.627722 1.438636 # [3,] 1.476511 1.627722
# 1.854880 1.620617 # [4,] 1.328423 1.438636 1.620617 1.984064 } else { vy_full
# <- NULL # ee_mean <- colMeans(ee,na.rm=T) # for later on individual i use to
# substract. # dont need this.  }

# # if (binary==T) v_of_mu=y_hat*(1-y_hat) else v_of_mu=rep(1,nrow(data)) #
# v_of_mu is the A in the formula, where Vi = Ai(1/2) Rwi Ai(1/2) v_of_mu <-
# switch(family, gaussian = rep(1,nrow(data)), binomial = y_hat*(1-y_hat),
# stop('Please specify a correct type of data distribution among 'gaussian' and
# 'binomial'!') ) data <- data.frame(data,X.intercept = 1,y_hat =
# y_hat,v_of_mu=v_of_mu,check.names=F) # we let data now include the col for
# intercept.


# ################################################ ##############construct the
# Rw################## ##############working correlation matrix########
# ################################################ #-------define the Rw for
# max_longiN, then reduced dimension can be subsetted out.--------------# # Rws
# <- foreach ( d_w = different_waves_we_have) %do% { # no need to parallel at
# all.  # Rws <- foreach ( d_w = max_longiN) %do% { # no need to parallel at all.
# # if (d_w > 1) { dim_R <- max_longiN # the nrow and ncol for Rw of this subject
# i.  if (dim_R > 1) { # only when in longitudinal setting if (corstr ==
# 'independence') { rho0 <- 0 } else rho0 <- fit0$geese$alpha # for cs and ar1, a
# scalar returned.  # for unstructured, a vector returned # e.g. Browse[2]>
# geefit_result$geese$alpha # alpha.1:2 alpha.1:3 alpha.1:4 alpha.2:3 alpha.2:4
# alpha.3:4 # 0.8511 0.7668 0.6997 0.8819 0.7828 0.8822


# #------building the Rw for subject i---------------# # work for no-missing data
# case.  independence_Rw = diag(1,nrow=dim_R,ncol=dim_R) ar1_Rw = rho0 ^
# abs(row(independence_Rw)-col(independence_Rw)) # for ar1 cs_Rw = rho0 ^ (
# abs(row(independence_Rw)-col(independence_Rw))>0 ) # for exchangeable #
# unstructured is a little bit more work unstructured_Rw = independence_Rw # need
# gdata package lowerTriangle(unstructured_Rw) <- rho0 for (j in
# 2:nrow(independence_Rw)) { for (i in 1:(j-1) ) { unstructured_Rw[i,j] <-
# unstructured_Rw[j,i] } } # now we have unstructured_Rw

# #-------give appropriate Rw for max_longiN------------------# Rw_full_dim =
# switch(corstr, independence = independence_Rw, ar1 = ar1_Rw, exchangeable =
# cs_Rw, unstructured = unstructured_Rw, stop('Please provide an explicit
# correlation structure among 'independence','exchangeable','ar1' or
# 'unstructured'!') ) } else Rw_full_dim <- NULL # } # else Rw = NULL #
# return(Rw) # } # names(Rws) <- max_longiN # each Rw within Rws has a name ==
# different_waves_we_have for subject i.  # Browse[4]> Rws # $`4` # [,1] [,2]
# [,3] [,4] # [1,] 1.0000000 0.8813135 0.7767136 0.6845282 # [2,] 0.8813135
# 1.0000000 0.8813135 0.7767136 # [3,] 0.7767136 0.8813135 1.0000000 0.8813135 #
# [4,] 0.6845282 0.7767136 0.8813135 1.0000000 Rws <- 1 # for cross-sectional
# data.  if (max_longiN > 1 & (!full_obs_flag)){

# Rws <- foreach ( msc = measurements_seq_categories ) %do% { # no need to
# parallel # if (nchar(msc) > 1) { dim_vec <- as.integer(
# strsplit(msc,split='')[[1]] ) dim_select_index <- full_dim_vec %in% dim_vec Rw
# = Rw_full_dim[dim_select_index,dim_select_index] # } else Rw = NULL return(Rw)
# } names(Rws) <- measurements_seq_categories # now we have the Rws storing all
# needed Rw types for data.  } #

# # ################################################ ################End:
# construct the Rw########### ################################################


# # #---- set up the empty container for U and I under null hypothesis.  #
# U<-matrix(0,ncol=1,nrow=ncol(x)) # the U full under null # #------AK is -EU(1
# derive); BK is EUU'.----------# # AK<-matrix(0,ncol=ncol(x),nrow=ncol(x)) #
# BK<-matrix(0,ncol=ncol(x),nrow=ncol(x)) # I <-
# matrix(0,ncol=ncol(x),nrow=ncol(x)) # # The cov for above U, same as BK should
# be. ncol robust for both df and matrix. length in matrix will be product of the
# nrow*ncol.


# ######################################### #------derive Uifor each nfam
# number----# # a func to derive the U vector for each subject i.  # this is for
# missingData part of the data.  derive_Ui <- function(data_chunk =
# subset(data,subset = eval(as.name(id)) == 1) ){ # browser() token_toGetRw <-
# paste(data_chunk[[waves ]], collapse='') v = data_chunk$v_of_mu # x_allPart =
# x[attr(x,'id') ==unique(data_chunk[[id]]),] # including intercept col, cov if
# any and par.  if (cov == 1) { # only intercept there.  x_allPart =
# data_chunk[c('X.intercept',par)] # including intercept col, cov if any and par.
# } else x_allPart = data_chunk[c('X.intercept',cov,par)] # including intercept
# col, cov if any and par.

# x_allPart <- as.matrix(x_allPart) # a matrix, e.g. 1-4 row.  dev <- x_allPart *
# v # the D in the U function # x_allPart * c(.5,.4,.6,.7) ==
# diag(c(.5,.4,.6,.7)) %*% x_allPart y_resid <- data_chunk$y_resid # residual
# vector, i.e. yi - yhati. length = k # if (exists('vys') ) { # we have vys in
# parent.frame iff commonVariance & family == 'gaussian' as specified in line.177
# # too slow for exists.  if (commonVariance & family == 'gaussian' ) { # we have
# vys in parent.frame iff commonVariance & family == 'gaussian' as specified in
# line.177 vy = vys[[token_toGetRw]] } else { # ee_mean_sub <-
# ee_mean[data_chunk[[waves ]]] # the corresponding ee_mean_sub part according to
# waves. # vy = y_resid %*% t(y_resid) # if not commonVariance }

# # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame
# # if nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.  Rw
# = Rws[[token_toGetRw ]] # Rws comming from parent.frame # use token_toGetRw to
# get the Rws corresponding element.

# # ################################################ # ##############construct
# the Rw################## commented, because this part of building Rw has
# migrated out of function, so did it just once for different dimension of
# subject i, e.g. 1~4.  # ##############working correlation matrix######## #
# ################################################ # #-------define the Rw if dim
# > 1.--------------# # if (nrow(data_chunk) > 1) { # dim_R <- nrow(data_chunk) #
# the nrow and ncol for Rw of this subject i.  # if (corstr == 'independence') {
# # rho0 <- 0 # } else rho0 <- fit0$geese$alpha # # for cs and ar1, a scalar
# returned.  # # for unstructured, a vector returned # # e.g. Browse[2]>
# geefit_result$geese$alpha # # alpha.1:2 alpha.1:3 alpha.1:4 alpha.2:3 alpha.2:4
# alpha.3:4 # # 0.8511 0.7668 0.6997 0.8819 0.7828 0.8822


# # #------building the Rw for subject i---------------# # # work for no-missing
# data case.  # independence_Rw = diag(1,nrow=dim_R,ncol=dim_R) # ar1_Rw = rho0 ^
# abs(row(independence_Rw)-col(independence_Rw)) # for ar1 # cs_Rw = rho0 ^ (
# abs(row(independence_Rw)-col(independence_Rw))>0 ) # for exchangeable # #
# unstructured is a little bit more work # unstructured_Rw = independence_Rw #
# need gdata package # lowerTriangle(unstructured_Rw) <- rho0 # for (j in
# 2:nrow(independence_Rw)) { # for (i in 1:(j-1) ) { # unstructured_Rw[i,j] <-
# unstructured_Rw[j,i] # } # } # # now we have unstructured_Rw

# # #-------give appropriate Rw for subject i------------------# # Rw =
# switch(corstr, # independence = independence_Rw, # ar1 = ar1_Rw, # exchangeable
# = cs_Rw, # unstructured = unstructured_Rw, # stop('Please provide an explicit
# correlation structure among 'independence','exchangeable','ar1' or
# 'unstructured'!') # ) # } # # #
# ################################################ # ################End:
# construct the Rw########### # ################################################

# #---------------start calculate U and Sigma of U------------------------# if
# (nrow(data_chunk) > 1) { # equivalent to nchar(token_toGetRw) > 1
# V<-diag(sqrt(v) )%*%Rw%*%diag(sqrt(v) ) # the V_i
# U0<-t(dev)%*%ginv(V)%*%y_resid # for one subject, U full function under the
# null.  #------AK is -EU(1 derive); BK is EUU'.----------# # if Rw() is
# correctly specified, then AK = BK AKi = t(dev) %*% ginv(V) %*% dev I0<-t(dev)
# %*% ginv(V) %*% vy %*% t(ginv(V)) %*% dev # assume common variance vy across
# samples. The Variance of above U0.  # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy)
# %*%t(ginv(V))%*%dev # not assume common variance across samples BKi = I0

# } else if (nrow(data_chunk)==1){ # all goes in to scalar.  V<-v^2 # scalar,
# since one trait, no R is needed.  U0<-matrix(dev,ncol=1)*c(y_resid)/V
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V BKi <- I0 AKi <-
# matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V } if (!commonVariance ) vy <- NULL #
# setting back to NULL for safty.  return(list(U = U0, I = I0, AK = AKi, BK =
# BKi) )


# } ######################################### #------derive Uifor each nfam
# number----# # a func to derive the U vector for each subject i.  # this is for
# complete part of the data.  derive_Ui_forFullObs <- function(data_chunk =
# subset(data,subset = eval(as.name(id)) == 1) ){ # browser() # token_toGetRw <-
# paste(data_chunk[[waves ]], collapse='') v = data_chunk$v_of_mu if (cov == 1) {
# # only intercept there.  x_allPart = data_chunk[c('X.intercept',par)] #
# including intercept col, cov if any and par.  } else x_allPart =
# data_chunk[c('X.intercept',cov,par)] # including intercept col, cov if any and
# par.  x_allPart <- as.matrix(x_allPart) # a matrix, e.g. 1-4 row.  dev <-
# x_allPart * v # the D in the U function # x_allPart * c(.5,.4,.6,.7) ==
# diag(c(.5,.4,.6,.7)) %*% x_allPart y_resid <- data_chunk$y_resid # residual
# vector, i.e. yi - yhati. length = k # if (exists('vys') ) { # we have vys in
# parent.frame iff commonVariance & family == 'gaussian' as specified in line.177
# # too slow for exists.  if (commonVariance & family == 'gaussian' ) { # we have
# vys in parent.frame iff commonVariance & family == 'gaussian' as specified in
# line.177 vy = vy_full } else { # ee_mean_sub <- ee_mean[data_chunk[[waves ]]] #
# the corresponding ee_mean_sub part according to waves. # vy = y_resid %*%
# t(y_resid) # if not commonVariance }

# # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame
# # if nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.  Rw
# = Rw_full_dim # Rws comming from parent.frame # use token_toGetRw to get the
# Rws corresponding element.  #---------------start calculate U and Sigma of
# U------------------------# if (nrow(data_chunk) > 1) { # equivalent to
# nchar(token_toGetRw) > 1 V<-diag(sqrt(v) )%*%Rw%*%diag(sqrt(v) ) # the V_i
# U0<-t(dev)%*%ginv(V)%*%y_resid # for one subject, U full function under the
# null.  #------AK is -EU(1 derive); BK is EUU'.----------# # if Rw() is
# correctly specified, then AK = BK AKi = t(dev) %*% ginv(V) %*% dev I0<-t(dev)
# %*% ginv(V) %*% vy %*% t(ginv(V)) %*% dev # assume common variance vy across
# samples. The Variance of above U0.  # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy)
# %*%t(ginv(V))%*%dev # not assume common variance across samples BKi = I0

# } else if (nrow(data_chunk)==1){ # all goes in to scalar.  V<-v^2 # scalar,
# since one trait, no R is needed.  U0<-matrix(dev,ncol=1)*c(y_resid)/V
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V BKi <- I0 AKi <-
# matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V } if (!commonVariance ) vy <- NULL #
# setting back to NULL for safty.  return(list(U = U0, I = I0, AK = AKi, BK =
# BKi) ) } # ########################################################### #
# #------derive Uifor each nfam number for full obs case----# # # (obsolete) a
# func to derive the U vector for each subject i when encounter full obs case.  #
# derive_Ui_forFullObs <- function(index = 1:nCluster) { # based on how many
# nCluster, i go through the length of v_of_mu...  # result <- foreach (i =
# index)%do%{ # individual number or family number.  # #
# muy1<-muy[(lag*(i-1)+1):(lag*i)] # 1-4, 5-8, 9-12. a vector #
# v=v_of_mu[(max_longiN*(i-1)+1):(max_longiN*i)] # 1-4, 5-8, 9-12. a vector for
# Ai # mux1<-as.matrix(x[(max_longiN*(i-1)+1):(max_longiN*i),]) # a matrix, e.g.
# 1-4 row extracted, with cols are cov and par. usually dim = c(k, 30) #
# dev<-mux1*v # the D in the U function #
# yy<-e[(max_longiN*(i-1)+1):(max_longiN*i)] # residual vector, i.e. yi - yhati.
# length = k # if (!commonVariance ) vy_full = yy %*% t(yy) # if not
# commonVariance # if (max_longiN>1){ # the trait # #
# V<-diag(sqrt(c(v)))%*%R%*%diag(sqrt(c(v))) # V<-diag(sqrt(v)
# )%*%Rw_full_dim%*%diag(sqrt(v) ) # given full data, only one Rw in Rws list.  #
# U0<-t(dev)%*%ginv(V)%*%yy # for one subject, U full function under the null.  #
# #------AK is -EU(1 derive); BK is EUU'.----------# # # if Rw() is correctly
# specified, then AK = BK # AKi = t(dev) %*% ginv(V) %*% dev

# # I0<-t(dev) %*% ginv(V) %*% vy_full %*% t(ginv(V)) %*% dev # assume common
# variance vy across samples. The Variance of above U0.  # #
# I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy) %*%t(ginv(V))%*%dev # not assume common
# variance across samples # BKi = I0 # }

# # if (max_longiN==1){ # all goes in to scalar.  # V<-v^2 # scalar, since one
# trait, no R is needed.  # U0<-matrix(dev,ncol=1)*c(yy)/V #
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy_full)/V/V # BKi <- I0 # AKi <-
# matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V # } # if (!commonVariance ) vy_full
# <- NULL # setting back to NULL for safty.

# # # U<-U+U0 # # I<-I+I0 # # AK = AK + AKi # # BK = BK + BKi # return(list(U =
# U0, I = I0, AK = AKi, BK = BKi) )

# # } # return(result)

# # }

# # # derive_Ui_cf <- cmpfun(derive_Ui) # not of much usage.  # #
# derive_Ui_forFullObs_cf <- cmpfun(derive_Ui_forFullObs)

# # ######################################################################## #
# ###########the below is original loop structure######################### #
# ###########confirmed the above func lead to same result as this######### #
# ######################################################################## #
# #------derive Uifor each nfam number----# # the for loop equivalent to above
# derive_Ui if (full_obs_flag) { # this below for loop will be faster than dlply
# operation on subsetting data.frame.  # U_results_all_subjects_collected <-
# derive_Ui_forFullObs(index = 1:nCluster) # apply on full index if (id == 'id')
# { U_results_all_subjects_collected <- dlply(.data = data, .variables =.(id),
# .fun =derive_Ui_forFullObs, .progress = 'text', .parallel = .parallel, .paropts
# = list(.packages=c('MASS','gdata','plyr')), .inform=F) } else
# U_results_all_subjects_collected <- dlply(.data = data, .variables
# =.(eval(as.name(id))), .fun =derive_Ui_forFullObs, .progress = 'text',
# .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F)


# # ###################################################### #
# ###################################################### #
# ###################################################### #
# ###################################################### } else { #----1st, we
# need to build the vector with NA in for nCluster.-----# # #
# full_vector_nClusterLong <- rep(1:4,times=nCluster) # insertNA <- fuction(x =
# waves_by_id[[1]] ){ # data[[waves]][] # full_dim_vec # } # as we dont need to
# put the whole data into dlply computations.  id_have_missingObs_index <-
# which(waves_by_id != max_longiN) # the index for id have missing measurements
# data_full_subset <- subset(data[],subset = !eval(id) %in%
# id_have_missingObs_index) # we dont need this in derive_Ui_forFullObs
# data_missingObs_subset <- subset(data[],subset = eval(id) %in%
# id_have_missingObs_index) # we do need this in derive_Ui # browser() if (id ==
# 'id') { U_results_full_subjects_collected_subset <- dlply(.data =
# data_full_subset, .variables =.(id), .fun =derive_Ui_forFullObs, .progress =
# 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')), .inform=F) } else
# U_results_full_subjects_collected_subset <- dlply(.data = data_full_subset,
# .variables =.(eval(as.name(id))), .fun =derive_Ui_forFullObs, .progress =
# 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F)

# # need to put before below v_of_mu, x, e reassign value.

# #------we need to re-build v_of_mu, x, e (which now exclude id chunks which has
# missing obs)--------------------# # v_of_mu <-
# data_missingObs_subset[['v_of_mu']] # e <- data_missingObs_subset[['y_resid']]
# # x <- data.frame(intercept = 1, subset(data_missingObs_subset, select =
# grep('X\\.',colnames(data_missingObs_subset),perl=T,value=T)), check.names=F)
# # in the same format as original x.  # attr(x,'id') <-
# data_missingObs_subset[eval(id)] # the corresponding id for each row in x #
# browser() # suppressWarnings( # Progress disabled when using parallel plyr #
# user system elapsed # 1.657 0.084 3.173 # Warning messages: # 1: <anonymous>:
# ... may be used in an incorrect context: .fun(piece, ...)

# # 2: <anonymous>: ... may be used in an incorrect context: .fun(piece, ...)  #
# ---confirmed, the above info is not any wrong in the code and process.----# #
# so use suppressWarnings to get rid of this info.  if (id == 'id') {
# U_results_missingObs_subjects_collected_subset <- dlply(.data =
# data_missingObs_subset, .variables =.(id), .fun =derive_Ui, .progress = 'text',
# .parallel = .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')),
# .inform=F) # U_results_all_subjects_collected <- dlply(.data = data, .variables
# =.(id), .fun =derive_Ui, .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')), .inform=F) } else {
# U_results_missingObs_subjects_collected_subset <- dlply(.data =
# data_missingObs_subset, .variables =.(eval(as.name(id))), .fun =derive_Ui,
# .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F) #
# U_results_all_subjects_collected <- dlply(.data = data, .variables
# =.(eval(as.name(id))), .fun =derive_Ui, .progress = 'text', .parallel =
# .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')),.inform=F) } # )
# U_results_all_subjects_collected <-
# c(U_results_full_subjects_collected_subset,U_results_missingObs_subjects_collected_subset)

# } # Browse[1]> names(U_results_all_subjects_collected[[1]]) # [1] 'U' 'I' 'AK'
# 'BK' U <- colSums(laply(U_results_all_subjects_collected,function(x) x$U ),
# na.rm=T) # na.rm =T is robust for missing values in par. such individual cell
# of missing value will be skipped from summation. (cellwise skip, so if only 1
# obs left, the colSums will be this 1 obs value.)  I <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$I ), na.rm=T) #
# verified the dimension is correct to sum over.  AK <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$AK ), na.rm=T)#
# verified the dimension is correct to sum over.  BK <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$BK ), na.rm=T)#
# verified the dimension is correct to sum over.


# #-------------# # alpha or beta, or gammma or beta. I.e., one for X (cov) one
# for Z (parameter).  ## calculate covariance matrix of gamma conditional on beta
# ## I11 is the covariance matrix of beta; I22 is the covariance matrix of
# gamma;I12 is the covariance matrix of beta-gamma

# # giving the right names in case in scalar case, the name will not be given
# automatically.  if (cov == 1) { covName <- 'X.intercept' } else { covName <-
# c('X.intercept',cov) } U=as.data.frame(U) # to transfrom possible scalar to
# matrix.  rownames(U)=c(covName,par) I=as.data.frame(I)
# dimnames(I)=list(c(covName,par),c(covName,par) ) AK = as.data.frame(AK) BK =
# as.data.frame(BK) dimnames(AK) = dimnames(BK) <-
# list(c(covName,par),c(covName,par) )

# #----------# Upar<-U[par,] # U of paramters (of interests) part.
# I11<-as.matrix(I[covName,covName]) # it's also BK
# I12<-t(as.matrix(I[par,covName])) I21<-as.matrix(I[par,covName])
# I22<-as.matrix(I[par,par]) # if ( !allowMisspecified) {
# V<-I22-I21%*%ginv(I11)%*%I12 } else if (allowMisspecified) { AK11 =
# as.matrix(AK[covName,covName]) AK12 = t(as.matrix(AK[par,covName])) AK21 =
# as.matrix(AK[par,covName]) AK22 = as.matrix(AK[par,par]) # V <- I22 - AK21 %*%
# ginv(AK11) %*% I12 - I21 %*% ginv(AK11) %*% AK12 + AK21 %*% ginv(AK11) %*% I11
# %*% ginv(AK11) %*% AK12 }

# #------------execute the standard gee score test--------------# score_statistic
# <- t(Upar)%*%ginv(V)%*%Upar # 38.19595 for correctly specified # 38.19595 for
# allow misspecified df=length(par) pval<-1-pchisq(score_statistic,df) # the gee
# score test pvalue.

# return(list(U=Upar,Cov=V,geeScoreTest_out=data.frame(df=df,
# Score.statistics=score_statistic,pval=pval))) }




### nCluster: ncluster,or n subject, or n familiy depending your data type.  data:
### the long format data with all information in, i.e. including covariates,
### variants (with parameter of interests), response variables.  id: specify the
### name of your cluster id, e.g. 'id' or 'subject' in colnames y: the response,
### can be index or name vector cov: the covariates (such as age, gender, etc), can
### be index or name vector default: 1, with only intercept, so populational mean.
### par: the variants(of which the parameters we are testing not equal to 0), can
### be index or name vector waves: variable name specifying the ordering of
### repeated mesurements on the same unit.  Also used in connection with missing
### values.the can be index or name vector corstr = allowCovNA: if enabled, data
### are filtered with addition to NAs in cov part (delete rowwisely.); if disabled,
### an error msg will be given if there are NAs encountered in cov part.
### allowParNA: if enabled, geescore will return normally with those NAs excluded
### in computation; if disabled, an error msg will be given if there are NAs
### encountered in par part.  allowYNA: if enabled, geescore will deal with
### unbalanced data; if disabled, an error msg will be given if there are NAs
### encountered in y part.  usage: source this script in main script.  Note: 1,For
### .parallel = TRUE, 'foreach' backend should be registered in main script, e.g.
### doMC, doSNOW, etc.  for this script, suggest easy task to disable .parallel
### since communication/split will take more time.  2, after modification in U <-
### colSums(laply(U_results_all_subjects_collected,function(x) x$U ), na.rm=T) #
### na.rm =T is robust for missing values in par. such individuals missing will be
### skipped from summation. Thus, Robust to NA in par.  if in y, the fit and
### residual will be N - n, where n is the # of missing rows. So ee matrix need to
### be remodified.  if in cov, geese/geeglm will report error whatever na.action
### equates.


# 3, dont let id = 'id', y = 'y', cov = 'cov', par = 'par' , waves = 'waves' in
# colnames(data), make necessary changes to keep the argument values are
# different to argument names. Otherwise, now the program may make error.  4,
# data are better given sorted in 'id' and 'waves' (maybe different names in
# colnames(data)). Though the program will also sort the data in this way.  5,
# commonVariance = T only make effect when family = ''gaussian.

# geescore_foreach<-function(nCluster = NULL,data,id,y,cov =
# 1,par,waves,corstr=c('independence','exchangeable','ar1','unstructured'),family
# = c('binomial','gaussian'), allowCovNA = TRUE , allowParNA = TRUE, allowYNA =
# TRUE, biasCorrection=T, commonVariance=T,.parallel=FALSE){
# ######################### ### obsolete arguments: # allowMisspecified=T #
# doesn't matter any more.  #########################

# # In geepack, corstr allow: # a character string specifying the correlation
# structure.  The # following are permitted: 'independence', 'exchangeable', #
# 'ar1', 'unstructured', 'userdefined', and 'fixed' # browser()
# stopifnot(is.data.frame(data)) # we need data to be data.frame instead of
# matrix.  data <- tbl_dt(data) # transform to data.table for speed.  # data <-
# na.omit(data) # the geese/geeglm program doesn't allow NA in y or cov part
# either, so if there is NA, delete rowwisely.  # if in y, the fit and residual
# will be N - n, where n is the # of missing rows. So ee matrix need to be
# remodified.  # if in cov, geese/geeglm will report error whatever na.action
# equates.  # if in par, already solved.  if (!is.numeric(data[,id,with=F]) ) {
# firstWavesSeq <- data[,list(N=.N),by=id][,N]
# setnames(data,id,sprintf('%s_originalBackup',id)) # rename invalid id.  newID
# <- rep(1:length(firstWavesSeq),times = firstWavesSeq) # data[,eval(id):=newID]
# # assign new id.  data[[eval(id)]] <- newID # assign new id.  } if (id == 'id')
# { data <- plyr::arrange(data, id, eval(as.name(waves)) ) # plyr allow expr in
# ...  } else data <- plyr::arrange(data, eval(as.name(id)), eval(as.name(waves))
# )

# #-----------# if (allowYNA){ y_NA_index <- with(data, is.na(eval(as.name(y))) )
# # y is always a single col, as required in long format.  data <-
# data[!y_NA_index,] # only filger NAs in y_NA_index.  } else
# tryCatch(na.fail(data[y]), error = stop('NA in y(response) are not allowed!
# Please either enable 'allowYNA' or filter/impute 'y' part by
# yourself!'),warning = identity) #---# # if (cov[1] != 1) { # in multipleTraits
# scenarios, we dont need this.  if (allowCovNA) { # cov_NA_index <- with(data,
# is.na(eval(as.name(cov))) ) cov_NA_index <-
# apply(data[,cov,with=F],1,function(x) any(is.na(x)) ) # for data.table.  data
# <- data[!cov_NA_index,] } else tryCatch(na.fail(data[,cov,with=F]), error =
# stop('NA in cov(covariates) are not allowed! Please either enable 'allowCovNA'
# or filter/impute 'cov' part by yourself!'),warning = identity) # } #---# #
# par_NA_index <- with(data, is.na(eval(as.name(par))) ) # single par, which
# possibly not gonna happen all the time.  # par_NA_index <- with(data,
# is.na(eval(as.name(par))) ) par_NA_index <-
# apply(data[,par,with=F],1,function(x) any(is.na(x)) ) # for data.table.  if
# (allowParNA & sum(par_NA_index)!=0 ) { # we didn't filter any row due to Par
# NAs.  message('There are NAs in the 'par' part, with 'allowParNA = TRUE', we
# will keep going as usual! (results are still reliable given small amount of
# NAs, it's your responsibility to evaluate how many NAs there are in 'par' part.
# Impute your 'par' part (usually they are genetic variants) if there are too
# many missings to keep the accuracy of the program.)') } else if (!allowParNA &
# sum(par_NA_index)!=0) { stop('NA in par(parameters of interests, usually
# genetic variants) are not allowed! Please either enable 'allowParNA' or impute
# 'par' part by yourself!') } # if sum(par_NA_index)==0, allowParNA doesn't
# matter any longer.  # we just go ahead to below part.

# #==============================# family <- match.arg(family) corstr <-
# match.arg(corstr) if (is.null(nCluster)) nCluster <-
# length(unique(data[[eval(id)]])) # for data.frame data, how many nCluster we
# have after all data filtering.

# ####################### #---------------------# # this part not really used in
# the below code, but can be useful. So keep it here.  waves_by_id <-
# table(data[[eval(id)]]) if (all (waves_by_id == max(waves_by_id) )) {
# full_obs_flag <- TRUE # perfect, we have full obs data, which will facilitate
# the below computation.  } else full_obs_flag <- FALSE


# #---------------------# #######################


# ####################### #---------------------# # this is used below for
# building Rw.  different_waves_we_have <- sort(unique(waves_by_id)) # e.g. 4 or
# [1] 1 2 3 4 depending on full obs or missing obs occurs.
# #---------------------# #######################



# # if (cov == 1) { # # only intercept there.  # x <- data.frame(intercept = 1,
# data[,par],check.names=F) # } else x<- data.frame(intercept = 1,
# data[,c(cov,par)],check.names=F) # for the latter dimension usage for whole U .
# # attr(x,'id') <- data[eval(id)] # the corresponding id for each row in x

# #------------fit the model under null--------------#
# form<-as.formula(paste(y,'~ 0 +',paste(cov,collapse='+'),sep='') ) # y is
# 'char' # cov <- c(1,'X.2') # for debug purpose.  # if (binary==T) geefit_result
# <- geeglm( form, id = eval(id),waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family='binomial') else geefit_result <- geeglm( form, id =
# eval(id),waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family='gaussian') # which also contains the $geese for above.  #
# if (binary==T) geefit_result <- geeglm( form, id = eval( as.name(id) ),waves =
# eval( as.name(waves) ),data = data, corstr=corstr,family='binomial') else
# geefit_result <- geeglm( form, id = eval( as.name(id) ),waves = eval(
# as.name(waves) ),data = data, corstr=corstr,family='gaussian') # which also
# contains the $geese for above.  if (id == 'id') { geefit_result <- geeglm(
# form, id = id,waves = eval( as.name(waves) ),data = data,
# corstr=corstr,family=family) # which also contains the $geese for above.  }
# else geefit_result <- geeglm( form, id = eval(as.name(id)),waves = eval(
# as.name(waves) ),data = data, corstr=corstr,family=family) # which also
# contains the $geese for above.  fit0 <- geefit_result# that is the model
# fitting under NULL hypothesis.  # scall <- match.call() # mnames <- c('', '',
# 'formula', 'data', 'offset', 'weights', 'subset', # 'na.action', 'id', 'waves',
# 'corp') # cnames <- names(scall) # cnames <- cnames[match(mnames, cnames, 0)] #
# mcall <- scall[cnames] # if (is.null(mcall$id)) # mcall$id <- as.name('id') #
# mcall[[1]] <- as.name('model.frame') # mcall[[2]] <- form # m <- eval(mcall,
# parent.frame()) # id <- model.extract(m, id) # waves <- model.extract(m,
# 'waves') #--------extract various infomation from model fit---------# #
# R=fit0$working.correlation # max_longiN <- max(data[eval(waves)] ) # the max
# measurements within one subject.  # max_longiN <-
# data[,max(eval(as.name(waves)))] # wrong, not correct for univariate trait
# data.  max_longiN <- max(waves_by_id) # the max measurements within one
# subject. # syntax for data.table obj.  full_dim_vec <- 1:max_longiN # for later
# usage in building Rws and vys.

# y_hat<-fitted(fit0) # only for binomial case usage, length = n x k, e.g. 36720
# # U1<-matrix(0,ncol=1,nrow=ncol(x)) # I<-matrix(0,ncol=ncol(x),nrow=ncol(x))
# e=fit0$residual # length = n x k, a.k.a the y_resid data <- mutate(data,y_resid
# = e)


# ## First get the whole score vector and information matrix under null if
# (full_obs_flag) { ee<-matrix(e,ncol=max_longiN,byrow=T) # residual data frame
# with dim = n * ni # here only for below vy.  # only work for non-missing data
# case.  } else { # ---a sub func here--- # # construct the ee residual matrix to
# n*ni.  make_ee_row <- function(data_chunk = subset(data,subset =
# eval(as.name(id)) == 1) ){ # browser() measurements_seq <- data_chunk[[waves]]
# y_resid_vector <- rep(NA, max_longiN) y_resid_vector[measurements_seq] <-
# data_chunk$y_resid return(y_resid_vector) } if (id == 'id') { ee <- daply(.data
# = data, .variables =.(id), .fun =make_ee_row, .progress = 'text', .parallel =
# .parallel, .paropts = list(.packages=c('MASS','gdata','plyr')), .inform=F) }
# else { ee <- daply(.data = data, .variables =.(eval(as.name(id))), .fun
# =make_ee_row, .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F) }

# } #===================================# # ---a sub func here--- # # collect all
# possible measurements_seq categories in the data.  # e.g. 1234, 123, 134,12,24,
# etc.  if (!full_obs_flag){ # only do it when we have missing data.
# collect_all_measurements_seqs <- function(data_chunk = subset(data,subset =
# eval(as.name(id)) == 1) ){ # browser() measurements_seq <- data_chunk[[waves]]
# return(measurements_seq) } if (id == 'id') { measurements_seq_pool <-
# dlply(.data = data, .variables =.(id), .fun =collect_all_measurements_seqs,
# .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')), .inform=F) } else {
# measurements_seq_pool <- dlply(.data = data, .variables =.(eval(as.name(id))),
# .fun =collect_all_measurements_seqs, .progress = 'text', .parallel = .parallel,
# .paropts = list(.packages=c('MASS','gdata','plyr')),.inform=F) }

# # measurements_seq_pool # will be used later in building subject specific Rw.
# measurements_seq_pool_pasted <- sapply(measurements_seq_pool,paste,collapse='')
# # for N subjects, length = n measurements_seq_categories <- names(
# table(measurements_seq_pool_pasted) ) }

# if (commonVariance & family == 'gaussian') { # # v-cov y is exactly v-cov e #
# vy_full<-cov(ee,use='pairwise.complete.obs') # k by k v-cov matrix;
# pairwise.complete.obs retains the most informations. k = max_longiN # vy_full<-
# matrix( rowSums(apply(ee,1,tcrossprod),na.rm=T)/ (nCluster),nrow =
# max_longiN,byrow=T) # k by k v-cov matrix; commonVariance only for quantitative
# case, binomial cannot adopt this. # biased, see this example show: # >
# apply(aaa,1,tcrossprod) # [,1] [,2] [,3] # [1,] 1 4 9 # [2,] 4 NA 18 # [3,] 7
# 16 27 # [4,] 4 NA 18 # [5,] 16 NA 36 # [6,] 28 NA 54 # [7,] 7 16 27 # [8,] 28
# NA 54 # [9,] 49 64 81 # > rowSums(apply(aaa,1,tcrossprod)) # [1] 14 NA 50 NA NA
# NA 50 NA 194 # > rowMeans(apply(aaa,1,tcrossprod)) # [1] 4.666667 NA 16.666667
# NA NA NA 16.666667 # [8] NA 64.666667 # > rowMeans(apply(aaa,1,tcrossprod),
# na.rm=T) # [1] 4.666667 11.000000 16.666667 11.000000 26.000000 41.000000
# 16.666667 # [8] 41.000000 64.666667 # > rowSums(apply(aaa,1,tcrossprod),
# na.rm=T)/ 3 # [1] 4.666667 7.333333 16.666667 7.333333 17.333333 27.333333
# 16.666667 # [8] 27.333333 64.666667 if (max_longiN == 1){ # or cross-sectional
# data.  vy_full <- mean(ee^2,na.rm=T) # scalar # Browse[2]>
# identical(as.vector(ee^2), apply(ee,1,tcrossprod)) # [1] TRUE

# } else vy_full<- matrix( rowMeans(apply(ee,1,tcrossprod),na.rm=T),nrow =
# max_longiN,byrow=T) # k by k v-cov matrix; commonVariance only for quantitative
# case, binomial cannot adopt this.  #---# # get various vy_partial for missing
# obs scenarios.  if (!full_obs_flag){ if (max_longiN == 1){ # or cross-sectional
# data.  vys <- vy_full } else { vys <- foreach ( msc =
# measurements_seq_categories ) %do% { # no need to parallel # if (nchar(msc) >
# 1) { # dont need this, as even one obs in either longitudinal missing data
# setting or cross-sectional setting, it has a vy comming from common vys.
# dim_vec <- as.integer( strsplit(msc,split='')[[1]] ) dim_select_index <-
# full_dim_vec %in% dim_vec vy = vy_full[dim_select_index,dim_select_index] # }
# else vy = NULL return(vy) } names(vys) <- measurements_seq_categories # now we
# have the Rws storing all needed Rw types for data.  } } # this cov func is the
# same as below. I.e., Sum_i { (ee[i, ] - ee_mean) %*% t(ee[i, ] - ee_mean) } /
# (N - 1). N is the nCluster and also the nrow(ee).  # ee_mean <- colMeans(ee) #
# ee_centered <- sweep(ee,2,ee_mean) # matrix(
# rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/ (nCluster-1 ),nrow =
# max_longiN,byrow=T) # the estimates of cov(r), r is the matrix of ri in a row.
# # Browse[2]> matrix( rowSums(apply(ee_centered,1,tcrossprod),na.rm=T)/
# (nCluster-1 ),nrow = max_longiN,byrow=T) # [,1] [,2] [,3] [,4] # [1,] 2.000849
# 1.731253 1.476511 1.328423 # [2,] 1.731253 1.978508 1.627722 1.438636 # [3,]
# 1.476511 1.627722 1.854880 1.620617 # [4,] 1.328423 1.438636 1.620617 1.984064

# # Browse[2]> vy # [,1] [,2] [,3] [,4] # [1,] 2.000849 1.731253 1.476511
# 1.328423 # [2,] 1.731253 1.978508 1.627722 1.438636 # [3,] 1.476511 1.627722
# 1.854880 1.620617 # [4,] 1.328423 1.438636 1.620617 1.984064 } else { vy_full
# <- NULL # ee_mean <- colMeans(ee,na.rm=T) # for later on individual i use to
# substract. # dont need this.  }

# # if (binary==T) v_of_mu=y_hat*(1-y_hat) else v_of_mu=rep(1,nrow(data)) #
# v_of_mu is the A in the formula, where Vi = Ai(1/2) Rwi Ai(1/2) v_of_mu <-
# switch(family, gaussian = rep(1,nrow(data)), binomial = y_hat*(1-y_hat),
# stop('Please specify a correct type of data distribution among 'gaussian' and
# 'binomial'!') ) # data <- mutate(data,X.intercept = 1,y_hat =
# y_hat,v_of_mu=v_of_mu) # we let data now include the col for intercept.  data
# <- mutate(data,y_hat = y_hat,v_of_mu=v_of_mu) # we dont need intercept, since
# they are already included in each row for each trait. (heterogeneous
# intercepts)


# ################################################ ##############construct the
# Rw################## ##############working correlation matrix########
# ################################################ #-------define the Rw for
# max_longiN, then reduced dimension can be subsetted out.--------------# # Rws
# <- foreach ( d_w = different_waves_we_have) %do% { # no need to parallel at
# all.  # Rws <- foreach ( d_w = max_longiN) %do% { # no need to parallel at all.
# # if (d_w > 1) { dim_R <- max_longiN # the nrow and ncol for Rw of this subject
# i.  if (dim_R > 1) { # only when in longitudinal setting if (corstr ==
# 'independence') { rho0 <- 0 } else rho0 <- fit0$geese$alpha # for cs and ar1, a
# scalar returned.  # for unstructured, a vector returned # e.g. Browse[2]>
# geefit_result$geese$alpha # alpha.1:2 alpha.1:3 alpha.1:4 alpha.2:3 alpha.2:4
# alpha.3:4 # 0.8511 0.7668 0.6997 0.8819 0.7828 0.8822


# #------building the Rw for subject i---------------# # work for no-missing data
# case.  independence_Rw = diag(1,nrow=dim_R,ncol=dim_R) ar1_Rw = rho0 ^
# abs(row(independence_Rw)-col(independence_Rw)) # for ar1 cs_Rw = rho0 ^ (
# abs(row(independence_Rw)-col(independence_Rw))>0 ) # for exchangeable #
# unstructured is a little bit more work unstructured_Rw = independence_Rw # need
# gdata package lowerTriangle(unstructured_Rw) <- rho0 for (j in
# 2:nrow(independence_Rw)) { for (i in 1:(j-1) ) { unstructured_Rw[i,j] <-
# unstructured_Rw[j,i] } } # now we have unstructured_Rw

# #-------give appropriate Rw for max_longiN------------------# Rw_full_dim =
# switch(corstr, independence = independence_Rw, ar1 = ar1_Rw, exchangeable =
# cs_Rw, unstructured = unstructured_Rw, stop('Please provide an explicit
# correlation structure among 'independence','exchangeable','ar1' or
# 'unstructured'!') ) } else Rw_full_dim <- NULL # } # else Rw = NULL #
# return(Rw) # } # names(Rws) <- max_longiN # each Rw within Rws has a name ==
# different_waves_we_have for subject i.  # Browse[4]> Rws # $`4` # [,1] [,2]
# [,3] [,4] # [1,] 1.0000000 0.8813135 0.7767136 0.6845282 # [2,] 0.8813135
# 1.0000000 0.8813135 0.7767136 # [3,] 0.7767136 0.8813135 1.0000000 0.8813135 #
# [4,] 0.6845282 0.7767136 0.8813135 1.0000000 Rws <- 1 # for cross-sectional
# data.  if (max_longiN > 1 & (!full_obs_flag)){

# Rws <- foreach ( msc = measurements_seq_categories ) %do% { # no need to
# parallel # if (nchar(msc) > 1) { dim_vec <- as.integer(
# strsplit(msc,split='')[[1]] ) dim_select_index <- full_dim_vec %in% dim_vec Rw
# = Rw_full_dim[dim_select_index,dim_select_index] # } else Rw = NULL return(Rw)
# } names(Rws) <- measurements_seq_categories # now we have the Rws storing all
# needed Rw types for data.  } #

# # ################################################ ################End:
# construct the Rw########### ################################################


# # #---- set up the empty container for U and I under null hypothesis.  #
# U<-matrix(0,ncol=1,nrow=ncol(x)) # the U full under null # #------AK is -EU(1
# derive); BK is EUU'.----------# # AK<-matrix(0,ncol=ncol(x),nrow=ncol(x)) #
# BK<-matrix(0,ncol=ncol(x),nrow=ncol(x)) # I <-
# matrix(0,ncol=ncol(x),nrow=ncol(x)) # # The cov for above U, same as BK should
# be. ncol robust for both df and matrix. length in matrix will be product of the
# nrow*ncol.


# ######################################### #------derive Uifor each nfam
# number----# # a func to derive the U vector for each subject i.  # this is for
# missingData part of the data.  derive_Ui <- function(data_chunk ){ # browser()
# token_toGetRw <- paste(data_chunk[[waves ]], collapse='') v =
# data_chunk$v_of_mu # x_allPart = x[attr(x,'id') ==unique(data_chunk[[id]]),] #
# including intercept col, cov if any and par.  # if (cov[1] == 1) { # # only
# intercept there.  # x_allPart = subset(data_chunk,select=c(par)) # including
# intercept col, cov if any and par. # dplyr syntax # # x_allPart =
# data_chunk[,c(par),with=F] # data.table syntax.  # } else # x_allPart =
# data_chunk[,c(cov,par),with=F] x_allPart = subset(data_chunk,select=c(cov,par))
# # in multipleTraits, we all need this, cov at least contains intercepts.


# x_allPart <- as.matrix(x_allPart) # a matrix, e.g. 1-4 row.  dev <- x_allPart *
# v # the D in the U function # x_allPart * c(.5,.4,.6,.7) ==
# diag(c(.5,.4,.6,.7)) %*% x_allPart y_resid <- data_chunk$y_resid # residual
# vector, i.e. yi - yhati. length = k # if (exists('vys') ) { # we have vys in
# parent.frame iff commonVariance & family == 'gaussian' as specified in line.177
# # too slow for exists.  if (commonVariance & family == 'gaussian' ) { # we have
# vys in parent.frame iff commonVariance & family == 'gaussian' as specified in
# line.177 vy = vys[[token_toGetRw]] } else { # ee_mean_sub <-
# ee_mean[data_chunk[[waves ]]] # the corresponding ee_mean_sub part according to
# waves. # vy = y_resid %*% t(y_resid) # if not commonVariance }

# # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame
# # if nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.  Rw
# = Rws[[token_toGetRw ]] # Rws comming from parent.frame # use token_toGetRw to
# get the Rws corresponding element.

# #---------------start calculate U and Sigma of U------------------------# if
# (nrow(data_chunk) > 1) { # equivalent to nchar(token_toGetRw) > 1
# V<-diag(sqrt(v) )%*%Rw%*%diag(sqrt(v) ) # the V_i
# U0<-t(dev)%*%ginv(V)%*%y_resid # for one subject, U full function under the
# null.  #------AK is -EU(1 derive); BK is EUU'.----------# ### AK and BK not
# needed any more, so comment to save time.  # if Rw() is correctly specified,
# then AK = BK # AKi = t(dev) %*% ginv(V) %*% dev I0<-t(dev) %*% ginv(V) %*% vy
# %*% t(ginv(V)) %*% dev # assume common variance vy across samples. The Variance
# of above U0.  # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy) %*%t(ginv(V))%*%dev # not
# assume common variance across samples # BKi = I0

# } else if (nrow(data_chunk)==1){ # all goes in to scalar.  V<-v^2 # scalar,
# since one trait, no R is needed.  U0<-matrix(dev,ncol=1)*c(y_resid)/V
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V # BKi <- I0 # AKi <-
# matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V } if (!commonVariance ) vy <- NULL #
# setting back to NULL for safty.  # return(list(U = U0, I = I0, AK = AKi, BK =
# BKi) ) # return(list(U = U0, I = I0 ) ) return(cbind(U = U0, I = I0 ) )


# } ######################################### #------derive Uifor each nfam
# number----# # a func to derive the U vector for each subject i. updated for
# using dt and dplyr # this is for complete part of the data.  #
# derive_Ui_forFullObs <- function(data_chunk = subset(data,subset =
# eval(as.name(id)) == 1) ){ # browser() derive_Ui_forFullObs <-
# function(data_chunk){ # browser() # token_toGetRw <- paste(data_chunk[[waves
# ]], collapse='') v = data_chunk$v_of_mu # v = data_chunk[,v_of_mu] # dt syntax
# # if (cov[1] == 1) { # # only intercept there.  # x_allPart =
# subset(data_chunk,select=c(par)) # including intercept col, cov if any and par.
# # dplyr syntax # # x_allPart = data_chunk[,c(par),with=F] # data.table syntax.
# # } else { # x_allPart = data_chunk[,c(cov,par),with=F] x_allPart =
# subset(data_chunk,select=c(cov,par)) # at least cov has intercept.1-4 # }
# x_allPart <- as.matrix(x_allPart) # a matrix, e.g. 1-4 row.  dev <- x_allPart *
# v # the D in the U function # x_allPart * c(.5,.4,.6,.7) ==
# diag(c(.5,.4,.6,.7)) %*% x_allPart y_resid <- data_chunk$y_resid # residual
# vector, i.e. yi - yhati. length = k # y_resid <- data_chunk[,y_resid] #
# residual vector, i.e. yi - yhati. length = k # dt syntax # if (exists('vys') )
# { # we have vys in parent.frame iff commonVariance & family == 'gaussian' as
# specified in line.177 # too slow for exists.  if (commonVariance & family ==
# 'gaussian' ) { # we have vys in parent.frame iff commonVariance & family ==
# 'gaussian' as specified in line.177 vy = vy_full } else { # ee_mean_sub <-
# ee_mean[data_chunk[[waves ]]] # the corresponding ee_mean_sub part according to
# waves. # vy = y_resid %*% t(y_resid) # if not commonVariance }

# # Rw = Rws[[as.character(nrow(data_chunk) ) ]] # Rws comming from parent.frame
# # if nrow ==1, Rw = NULL, below code will not use Rw at all, so works fine.  Rw
# = Rw_full_dim # Rws comming from parent.frame # use token_toGetRw to get the
# Rws corresponding element.  #---------------start calculate U and Sigma of
# U------------------------# if (nrow(data_chunk) > 1) { # equivalent to
# nchar(token_toGetRw) > 1 V<-diag(sqrt(v) )%*%Rw%*%diag(sqrt(v) ) # the V_i
# U0<-t(dev)%*%ginv(V)%*%y_resid # for one subject, U full function under the
# null.  #------AK is -EU(1 derive); BK is EUU'.----------# # if Rw() is
# correctly specified, then AK = BK # AKi = t(dev) %*% ginv(V) %*% dev I0<-t(dev)
# %*% ginv(V) %*% vy %*% t(ginv(V)) %*% dev # assume common variance vy across
# samples. The Variance of above U0.  # I0<-t(dev)%*%ginv(V)%*% yy %*% t(yy)
# %*%t(ginv(V))%*%dev # not assume common variance across samples # BKi = I0

# } else if (nrow(data_chunk)==1){ # all goes in to scalar.  V<-v^2 # scalar,
# since one trait, no R is needed.  U0<-matrix(dev,ncol=1)*c(y_resid)/V
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V # BKi <- I0 # AKi <-
# matrix(dev,ncol=1)%*%matrix(dev,nrow=1)/V } if (!commonVariance ) vy <- NULL #
# setting back to NULL for safty.  # return(list(U = U0, I = I0, AK = AKi, BK =
# BKi) ) # return(list(U = U0, I = I0) ) return(cbind(U = U0, I = I0) ) } #
# ########################################################### #---define a very
# useful function na+ here-------------------# `na+` <- function(x,y) {
# x[is.na(x)] <- 0; y[is.na(x)] <- 0; x + y # this allow pairwise +, not sum to a
# scalar.  }


# # ######################################################################## #
# ###########the below is original loop structure######################### #
# ###########confirmed the above func lead to same result as this######### #
# ######################################################################## #
# #------derive Uifor each nfam number----# # the for loop equivalent to above
# derive_Ui # data_by_id <- grouped_dt(data,quote(id)) # produce grouped_dt, for
# do() usage.  # browser() if (full_obs_flag) { # this below for loop will be
# faster than dlply operation on subsetting data.frame.  #
# U_results_all_subjects_collected <- derive_Ui_forFullObs(index = 1:nCluster) #
# apply on full index # setkeyv(data,id) # data_chunks <- iter(data,by=column)
# ids <- unique(data[[id]]) U_results_all_subjects_collected <- foreach(id_s =
# ids, .combine='na+') %dopar%{ derive_Ui_forFullObs(data_chunk =
# filter(data,eval(as.name(id)) %in% id_s) ) } # } # #
# U_results_all_subjects_collected <- do(data_by_id,derive_Ui_forFullObs) # do is
# slow, not fast at all in dplyr.  # } else # U_results_all_subjects_collected <-
# dlply(.data = data, .variables =.(eval(as.name(id))), .fun
# =derive_Ui_forFullObs, .progress = 'text', .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F) #
# U_results_all_subjects_collected <- dlply(.data = data, .variables
# =.(eval(as.name(id))), .fun =derive_Ui_forFullObs, .progress = 'text',
# .parallel = .parallel, .paropts =
# list(.packages=c('MASS','gdata','plyr')),.inform=F)


# # ###################################################### #
# ###################################################### #
# ###################################################### #
# ###################################################### } else { #----1st, we
# need to build the vector with NA in for nCluster.-----# # #
# full_vector_nClusterLong <- rep(1:4,times=nCluster) # insertNA <- fuction(x =
# waves_by_id[[1]] ){ # data[[waves]][] # full_dim_vec # } # as we dont need to
# put the whole data into dlply computations.

# id_have_missingObs_index <- as.integer( names(waves_by_id)[which(waves_by_id !=
# max_longiN)] )# the index for id have missing measurements # data_full_subset
# <- subset(data,subset = !eval(id) %in% id_have_missingObs_index)
# data_full_subset <- data[!eval(as.name(id))%in%id_have_missingObs_index,] #
# data_missingObs_subset <- subset(data,subset = eval(id) %in%
# id_have_missingObs_index) data_missingObs_subset <-
# data[eval(as.name(id))%in%id_have_missingObs_index,] browser()
# #---------------------------# # a, for data_full_subset ids <-
# unique(data_full_subset[[id]]) nWorkers <- getDoParWorkers() approxSingleN <-
# round(length(ids)/nWorkers) clusterExport(cl,'naPlus') system.time(
# U_results_full_subjects_collected_subset <- foreach(id_s = iter(t(ids[1:1500]),
# chunks=approxSingleN),
# .combine='naPlus',.verbose=F,.errorhandling='stop',.packages=c('geepack','MASS','doSNOW','plyr','dplyr','gdata'),.inorder=FALSE)
# %dopar%{ foreach(id_ss = id_s, .combine='naPlus',.inorder=FALSE) %do% {
# derive_Ui_forFullObs(data_chunk =
# filter(data_full_subset,data_full_subset[[id]] %in% id_ss) ) } # dim(
# filter(data_full_subset,eval(as.name(id)) %in% id_s) ) # doesn't work in
# parallel env, does work in serial.  # dim(
# filter(data_full_subset,data_full_subset[[id]] %in% id_s) ) # work } ) # 130s #
# compared with system.time( U_results_full_subjects_collected_subset <-
# dlply(.data = data_full_subset[1:(1500*4)], .variables =.(eval(as.name(id))),
# .fun =derive_Ui_forFullObs, .progress = 'text', .parallel = F, .paropts =
# list(.packages=c('MASS','gdata','plyr')), .inform=F) ) # 15s system.time(
# U_results_full_subjects_collected_subset <- by(data =
# data_full_subset[1:(1500*4)], INDICES = data_full_subset[1:(1500*4),id,with=F],
# FUN =derive_Ui_forFullObs) ) # 42.580 # system.time(
# U_results_full_subjects_collected_subset <- tapply(X =
# data_full_subset[1:(1500*4)],
# unlist(data_full_subset[c(1:(1500*4)),id,with=F]), derive_Ui_forFullObs) ) #
# 15s

# # U<-matrix(0,ncol=1,nrow=ncol(x)) I<-matrix(0,ncol=length(x),nrow=length(x))
# d_f_subset <- data_full_subset[1:(1500*4)] muy <- d_f_subset$y_hat v_of_mu <-
# d_f_subset$v_of_mu y_resid <- d_f_subset$y_resid lag=max_longiN x <-
# d_f_subset[,c(labels_x_new,labels_predictors_new),with=F] system.time( for (i
# in 1: (nrow(d_f_subset)/max_longiN)) { muy1<-muy[(lag*(i-1)+1):(lag*i)]
# v=v_of_mu[(lag*(i-1)+1):(lag*i)] mux1<-as.matrix(x[(lag*(i-1)+1):(lag*i),])
# dev<-mux1*v yy<-y_resid[(lag*(i-1)+1):(lag*i)]

# if (lag>1){ V<-diag(sqrt(c(v)))%*%Rw_full_dim%*%diag(sqrt(c(v)))
# U0<-t(dev)%*%ginv(V)%*%yy I0<-t(dev)%*%ginv(V)%*%vy_full%*%t(ginv(V))%*%dev }

# if (lag==1){ V<-v^2 U0<-matrix(dev,ncol=1)*c(yy)/V
# I0<-matrix(dev,ncol=1)%*%matrix(dev,nrow=1)*c(vy)/V/V }

# U<-U+U0 I<-I+I0 } ) #


# #---------------------------# # b, for data_missingObs_subset ids <-
# unique(data_missingObs_subset[[id]]) system.time(
# U_results_missingObs_subjects_collected_subset <- foreach(id_s = ids,
# .combine='na+') %dopar%{ derive_Ui(data_chunk =
# filter(data_missingObs_subset,eval(as.name(id)) %in% id_s) ) } )


# # need to put before below v_of_mu, x, e reassign value.

# #------we need to re-build v_of_mu, x, e (which now exclude id chunks which has
# missing obs)--------------------# # v_of_mu <-
# data_missingObs_subset[['v_of_mu']] # e <- data_missingObs_subset[['y_resid']]
# # x <- data.frame(intercept = 1, subset(data_missingObs_subset, select =
# grep('X\\.',colnames(data_missingObs_subset),perl=T,value=T)), check.names=F)
# # in the same format as original x.  # attr(x,'id') <-
# data_missingObs_subset[eval(id)] # the corresponding id for each row in x #
# browser() # suppressWarnings( # Progress disabled when using parallel plyr #
# user system elapsed # 1.657 0.084 3.173 # Warning messages: # 1: <anonymous>:
# ... may be used in an incorrect context: .fun(piece, ...)

# # 2: <anonymous>: ... may be used in an incorrect context: .fun(piece, ...)  #
# ---confirmed, the above info is not any wrong in the code and process.----# #
# so use suppressWarnings to get rid of this info.  # )
# U_results_all_subjects_collected <-
# `na+`(U_results_full_subjects_collected_subset,U_results_missingObs_subjects_collected_subset)

# } # Browse[2]> U_results_all_subjects_collected %>% dim # [1] 212 213

# # Browse[1]> names(U_results_all_subjects_collected[[1]]) # [1] 'U' 'I' 'AK'
# 'BK'

# ################################################## #------very
# slow-------------- # since large matrix manipulation in from ever I step
# calculation # U <- colSums(laply(U_results_all_subjects_collected,function(x)
# x$U ), na.rm=T) # na.rm =T is robust for missing values in par. such individual
# cell of missing value will be skipped from summation. (cellwise skip, so if
# only 1 obs left, the colSums will be this 1 obs value.)  # I <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$I ), na.rm=T) #
# verified the dimension is correct to sum over.  # AK <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$AK ), na.rm=T)#
# verified the dimension is correct to sum over.  # BK <-
# colSums(laply(U_results_all_subjects_collected,function(x) x$BK ), na.rm=T)#
# verified the dimension is correct to sum over.
# ################################################## U <-
# U_results_all_subjects_collected[,1] I <- U_results_all_subjects_collected[,-1]

# # #----use it on getting the U, I, AK and BK------------# # U <-
# lapply(U_results_all_subjects_collected,function(x) x$U) # U <- Reduce(`na+`,U)

# # I <- lapply(U_results_all_subjects_collected,function(x) x$I) # I <-
# Reduce(`na+`,I)


# #---------------------------------# ### since the result is the same whatever
# the allowMisspecified is, ### we will not use AK and BK any more.  # AK <-
# lapply(U_results_all_subjects_collected,function(x) x$AK) # AK <-
# Reduce(`na+`,AK)

# # BK <- lapply(U_results_all_subjects_collected,function(x) x$BK) # BK <-
# Reduce(`na+`,BK)


# # browser() #

# #-------------# # alpha or beta, or gammma or beta. I.e., one for X (cov) one
# for Z (parameter).  ## calculate covariance matrix of gamma conditional on beta
# ## I11 is the covariance matrix of beta; I22 is the covariance matrix of
# gamma;I12 is the covariance matrix of beta-gamma

# # giving the right names in case in scalar case, the name will not be given
# automatically.  # if (cov[1] == 1) { # covName <- 'X.intercept' # } else {
# covName <- cov # we only have one case, since intercept is already integrated
# in cov in multipleTraits scenarios # } U=as.matrix(U) # to transfrom possible
# scalar to matrix.  rownames(U)=c(covName,par) I=as.matrix(I)
# dimnames(I)=list(c(covName,par),c(covName,par) ) # AK = as.matrix(AK) # BK =
# as.matrix(BK) # dimnames(AK) = dimnames(BK) <-
# list(c(covName,par),c(covName,par) )

# #----------# Upar<-U[par,,drop=F] # U of paramters (of interests) part.
# I11<-as.matrix(I[covName,covName]) # it's also BK
# I12<-t(as.matrix(I[par,covName])) I21<-as.matrix(I[par,covName])
# I22<-as.matrix(I[par,par]) # # if ( !allowMisspecified) {
# V<-I22-I21%*%ginv(I11)%*%I12 # } else if (allowMisspecified) { # AK11 =
# as.matrix(AK[covName,covName]) # AK12 = t(as.matrix(AK[par,covName])) # AK21 =
# as.matrix(AK[par,covName]) # AK22 = as.matrix(AK[par,par]) # # # V <- I22 -
# AK21 %*% ginv(AK11) %*% I12 - I21 %*% ginv(AK11) %*% AK12 + AK21 %*% ginv(AK11)
# %*% I11 %*% ginv(AK11) %*% AK12 # }

# #------------execute the standard gee score test--------------# score_statistic
# <- t(Upar)%*%ginv(V)%*%Upar # 38.19595 for correctly specified # 38.19595 for
# allow misspecified df=length(par) pval<-1-pchisq(score_statistic,df) # the gee
# score test pvalue.

# return(list(U=Upar,Cov=V,geeScoreTest_out=data.frame(df=df,
# Score.statistics=score_statistic,pval=pval))) }
 
