#!/work/02040/yyang/Software/R/bin/Rscript --slave --vanilla
# Content: the main executable to run association test using aSPU methods.
# Author: Yang Yang
# Demo arguments:
# filename_genotype <- "genotpe_data_sample.Rdata"
# filename_annotation <- "genotype_annotation.Rdata"
# filename_phenotype <- "phe_cov_longitudinal_longFormat_sampleData.Rdata"
# responseVarName <- "trigs_longi"
# covariateVarFile <- "covariateNames.txt"
# workingCorstr_input <- "independence"
# b_times <- 1e3
# seed <- 123456
# usePermutedU <- FALSE
# parallel_scheme <- "SNOW"
# jobExecuteOnMultipleNodes <- FALSE
# prefix <- "singleNode_simulation"
# include_aSPU <- TRUE
# parallel_over_gene <- FALSE
###########################################################################
suppressPackageStartupMessages(require(getopt))	# load necessary packages
 spec = matrix(c(
   'verbose', 'v', 0, "logical","display verbose message printed",
   'help'   , 'h', 0, "logical","call 'help' document",
   'filename_genotype'  , 'f', 1, "character","the Rdata containing the genotype data you want to run program on as explanatory variables",
   'filename_annotation'  , 'a', 1, "character","the Rdata containing the genotype data annotation",
   'filename_phenotype'  , 'P', 1, "character","the Rdata containing the phenotype information and all covariate information as columns in a df/matrix.",
   # 'chr'  , 'c', 1, "integer","the chromosome # from 1~25, dont put X or Y please.",
   'responseVarName'  , 'r', 1, "character","specify the response variable name as the column in 'filename_phenotype' object, e.g. trigs_longi",
   'covariateVarFile'  , 'c', 2, "character","Optional. When you want to adjust for covariates in the association study, specify the text filename (e.g. covariateNames.txt), which stores covariate variable names in the format of one name per line. The covariate info are stored together with responseVarName as in columns in 'filename_phenotype' object. See example file. ",
   'workingCorstr_input'   , 'w', 1, "character","the Rw you want to use among 'independence', 'exchangeable', 'ar1' and 'unstructured'",
   'include_aSPU' , 'u', 1, "logical","TRUE: the analysis will execute default 5 tests and include aSPU family tests, which cost more time; FALSE: the analysis will only execute 5 tests (sum, SSU, SSUw, Score and UminP).",
	 'b_times'  , 'b', 1, "integer","the iterations for both simulation or permutation, recommend to start with 1e3 for first time genome-wide scan; then use 1e6 for only those significant genes in the initial scan stage (you only need to update the 'filename_annotation' to include only those significant genes). This option is only effective to aSPU family tests.",
   'seed'  , 's', 2, "integer","Optional. You can specify a random integer, say 1234, so that the result can be repeated in the future if you rerun the program on the same input files. The default seed is 123. This option is only effective to aSPU family tests",
   'usePermutedU' , 'U', 1, "logical","TRUE: use permutation method to get U null distribution; FALSE: use simulation method to get U null distribution according to MVN distribution. This option is only effective to aSPU family tests.",
   'parallel_scheme'     , 'S', 1, "character","SERIAL, SNOW or MC (the latter two are parallel schemes. SNOW is memory efficient while communications among nodes cost heavily; MC is communications efficient while taking up much more memory). Recommend to use 'SERIAL', when you have a small sample size and a few genes; use 'MC' when your memory is big enough (e.g. single node with >= 32gb mem); use 'SNOW' when node memory is limited.",
   'parallel_over_gene'     , 'g', 1, "logical","TRUE: parallel computing over genes; FALSE: parallel computing within each gene while over genes serial computing is adopted. For small samples (e.g. <= 3000 subjects), recommend to use 'TRUE'; for larger samples, recommend to use 'FALSE'. If 'usePermutedU' = TRUE, always set it to 'FALSE'. This option doesn't matter when you use 'SERIAL' as the 'parallel_scheme'. ",
   'jobExecuteOnMultipleNodes'     , 'j', 1, "logical","TRUE of FALSE (parallel run on single node or multiple nodes. Only 'SNOW' supports the latter.)",
   'MultipleNodes_profile','m',2,'character',"This options comes together with 'jobExecuteOnMultipleNodes'==TRUE. The file contains the configuration for the virtual cluster with multiple nodes you intent to use. Check the example file.",
   'prefix', 'e', 2, "character","Optional. It specifies the prefix you want to append to any output file names, e.g. specify 'Prefix1213' to make the output filename to be 'Prefix1213_Chromosome1ExomeChip.Rdata'"
 ), byrow=TRUE, ncol=5);
 opt = getopt(spec);
 print(opt)
 # if help was asked for print a friendly message
 # and exit with a non-zero error code
 if ( length(opt) == 1 | !is.null(opt$help) ) {	# conditions call 'help'
   cat(getopt(spec, usage=TRUE));
   q(status=1);
 }

 #print some progress messages to stderr, if requested.
 if ( !is.null(opt$verbose) ) { write("Parsing commands...",stdout()); }

 for(i in 1:length(opt))  {
	assign(names(opt)[i],opt[[i]])	# assign opt member obj into .GlobalEnv for later call usage.
 }

#---------record the absolute working root path, where the main executable locates.-----------# 
programRootPWD <-  getwd() 	# for later usage of root path.
#---------record the OutputData path, where all the results produced go there.-----------# 
outputPWD <- sprintf("%s/OutputData/",programRootPWD)
#
#---------setup the dumpto file for collecting errors---------------#
if ( exists("prefix") ) {
	options(error = quote({
	message("Please check log file under program root and error dump file in", outputPWD );
	dump.frames(dumpto = sprintf("%s/%s.last.dump",outputPWD,prefix),to.file = TRUE); 
	q()}))	# for debugging purpose, will output a file "$prefix.last.dump.rda" for further debugging if error occured.
} else 	{
	prefix <- "Default"	# give the default value.
	options(error = quote({
	message("Please check log file under program root and error dump file in", outputPWD );
	dump.frames(dumpto = sprintf("%s/%s.last.dump",outputPWD,prefix),to.file = TRUE); 
	q()}))	# for debugging purpose, will output a file "last.dump.rda" for further debugging if error occured.
}
#===================================================================================#
# start logging.
log_connection <- file(sprintf("%s/%s_aSPU_associationTest.log",programRootPWD,prefix), open = "wt")	# log file always put under './' deliberately (it's convenient after you run the program, you dont need to change folder to look at the log file.)
sink(log_connection, type="output",split = TRUE)
sink(log_connection, type="message",split = FALSE)	# cannot split, so prefer save to log for logging the errors.
# on.exit( catt(geterrmessage()) )	# only used in function.
# the example we use is
cat("=====================================================\n")

#---#
# transforming jobExecuteOnMultipleNodes to another variable of class "logical".
if ( grepl("^T",jobExecuteOnMultipleNodes,perl=T,ig=T) ) jEOM <- T else if ( grepl("^F",jobExecuteOnMultipleNodes,perl=T,ig=T) ) jEOM = F else {    
	cat(getopt(spec, usage=TRUE));
	cat('========================================================\n')
	stop("Something wrong with parallel_scheme appointment!")
    q(status=1);
}

cat("=====================================================\n")
cat("Start at ", date(),"\n" )
cat("parallel_scheme choice is:",parallel_scheme,"\n")
# the example we use is
cat("=====================================================\n")
suppressPackageStartupMessages(source("head.r") ) 
#----------------setting up the parallel env------------------------#
# ---a subfunc to setup the parallel env---#
setParallel <- function(parallel_scheme = parallel_scheme){
	switch(parallel_scheme,
			MC = registerDoMC(detectCores()),
			SNOW = {cl <<- makeSOCKcluster(detectCores())	# make it mutable object.
					setDefaultCluster(cl)
					registerDoSNOW(cl)
					#---need to send the path and point to the correct path on each slave node---
					clusterExport(cl, "programRootPWD")
					clusterEvalQ(cl,source( sprintf('%s/head.r',programRootPWD) )) # load necessary packages on each, slave must be given absolute path, REMEBER!
					# on.exit({stopCluster(cl)})	# turn if off savely, not put here, otherwise generate -> kill immediately.
					},
			SERIAL = registerDoSEQ() ,
			stop("Not valid parallel_scheme appointed!")
	)
	#
	catt("We have registered ",getDoParWorkers()," cpu(s) on machine: ",system("hostname",intern=T), " by ", shQuote(getDoParName()) )
}
#-------#
if (!jEOM) {
	setParallel(parallel_scheme)
	on.exit({stopCluster(cl)})
} else if (jEOM && exists("MultipleNodes_profile") ) {
	m_n_profile <- read.table(sprintf("%s/NodesProfile/%s",programRootPWD,MultipleNodes_profile),header=T)
	nodes_file <- unlist( apply(m_n_profile, 1, function(x) rep(x['nodeName'],x['cpus'])) )	# the cpus pool.
	cl <- makeSOCKcluster(nodes_file)
	setDefaultCluster(cl)
	registerDoSNOW(cl)
	#---need to send the path and point to the correct path on each slave node---
	clusterExport(cl, "programRootPWD")
	clusterEvalQ(cl,source( sprintf('%s/head.r',programRootPWD) )) # load necessary packages on each, slave must be given absolute path, REMEBER!
	on.exit({stopCluster(cl)})	# turn if off savely.
	catt("We have registered ",getDoParWorkers()," cpu(s) by ", shQuote(getDoParName()), "on machine: " )
	print(m_n_profile)
} else	{
	setParallel(parallel_scheme)	# still run on single node.
	on.exit({stopCluster(cl)})
}
# # clusterEvalQ(cl,ls())	# for debugging purpose.
catt("Now we have ",getDoParWorkers(),  " cores running!")

#====================source the kernal functions==========================#
source("transform.r")
source("aspuTest_parralled_withSco_RV_version_updatedScoreTest_within.r")	# load SPU3 func.
source("gee_scoreTest_self_geepack_allowMissingData_forMultiTraitsExt_needCovWithInterceptCol_trial.r") # load geescore_parallel func, very slow
source("RV_permutation_longi.r") # load RV_permute_test
source("Tests5V1.r")	# execute Dr. Pan Wei's Score, SSU, SSUw, Sum, UminP tests.
#---send all functions to all nodes if cl exists---#
if (exists("cl")) clusterEvalQ(cl, {
															setwd(programRootPWD)
															source("transform.r")
															source("aspuTest_parralled_withSco_RV_version_updatedScoreTest_within.r")	# load SPU3 func.
															source("gee_scoreTest_self_geepack_allowMissingData_forMultiTraitsExt_needCovWithInterceptCol_trial.r") # load geescore_parallel func, very slow
															source("RV_permutation_longi.r") # load RV_permute_test
															source("Tests5V1.r")	# execute Dr. Pan Wei's Score, SSU, SSUw, Sum, UminP tests.
															})
# # clusterEvalQ(cl,ls())	# for debugging purpose.
###############################################
#------load data from here--------------------#
setwd("GenotypeData")	# the folder storing the genotype data.
#---SNV part---#
real_Data <- load(filename_genotype)	# the real data input for chr x.
# GT
GT <- eval(as.name(real_Data))	# store the genotype data in GT variable
rownamesGT <- rownames(GT)
MAF_all <-apply(X=GT, 2, mean,na.rm=T)	# to QC on MAF across all variants.
setwd("..")
#------load the snp annotation file---#
setwd("AnnotationData")	# the folder storing the genotype data.
real_Data_annotation <- load(filename_annotation)	# constant annotation file.
snpinfo <- eval(as.name(real_Data_annotation))	# store the genotype data in GT variable
setwd("..")

#---Phenotype part and covariates part---#
setwd("PhenotypeCovariateData")
real_Data_phenotype <- load(filename_phenotype)	# the file contains long format for EA + AA.
phe_cov_longitudinal_longFormat <- eval( as.name( real_Data_phenotype))
# > phe_cov_longitudinal_longFormat %>% dim
# [1] 43572    18
 # phe_cov_longitudinal_longFormat$PCA_V11 %>>% (is.na(.)) %>>% sum
# [1] 72
setwd("..")

#=========start to join the data parts====================#
# ---transform to dt instead of df---
GT <- as.data.table(GT,keep.rownames=T)	# rownames as rn
# AIMs_annotation <- as.data.table(AIMs_annotation)
snpinfo <- as.data.table(snpinfo)
phe_cov_longitudinal_longFormat <- as.data.table(phe_cov_longitudinal_longFormat)

# ---keep those only appearing in annotation file 'snpinfo'
commonSNVNames <- intersect( colnames(GT),snpinfo$Name )
snpinfo <- subset(snpinfo, subset = Name %in% commonSNVNames)	# update the content of snpinfo
GT <- GT[, c("rn",commonSNVNames),with=FALSE]	# rn need to be in.

# ---set the key for phe_cov data object---#
setkeyv(phe_cov_longitudinal_longFormat,"id")


#############################################################################
#---subfunc: define a gene name based iter object and corresponding method---.
 iterByGene <- function(a, geneNames,annotation, key="single_gene",patientID = "rn") {	# need a and annotation be data.table, geneNames be character vector.
	setkeyv(annotation,key)	# set key
	i <- 1
  nextElem <- function() {
	geneName <- geneNames[i]	#i will change
    if ( is.na(geneName) ) stop('StopIteration')
    SNV_names <- annotation[geneName,Name]
    i <<- i + 1
	a[,c(patientID,SNV_names), with = FALSE]
  }

  structure( list(nextElem=nextElem), class=c('iblkcol', 'iter') )	# so we will have the S3 method "nextElem.iblkcol" applicable, which is loaded in header file.
}

#---execute---
allGeneNames <- snpinfo$single_gene %>>% unique %>>% na.omit
iterObj <- iterByGene(a = GT, geneNames = allGeneNames, annotation = snpinfo, key = "single_gene" )

#----get the covNames ready, the Z part in H*theta-----#
if (exists("covariateVarFile")) {
	covariateVar <- fread(sprintf( "%s/CovariateNames/%s", programRootPWD, covariateVarFile ), header=F)
	covNames <- c("intercept",unlist(covariateVar) )	# 1st two is enough as from observation on file "white.eval".
} else {
	covNames <- c("intercept")
}


if (!include_aSPU) {	# only 5 asymptotic tests, i.e., sum, SSU, SSUw, Score and UminP.
		parallelFlag <- FALSE
			#--run---#
		testResult <- foreach(iter = iterObj,.errorhandling = "pass",.verbose=T) %dopar% {	# since asymptotic tests are fast, we prefer to parallel over genes.
		 # registerDoMC(detectCores())	# for each in cl_mother
		 # iter = nextElem(iterObj)	# for debugging usage.
		 #----some limited QC here, more extensive QC should be done outside of this program----#
		 # ---1, NA control----#
		 NaRatios <- apply(iter,2,function(x) mean(is.na(x)) )
		 iter <- iter[,NaRatios <= 0.05,with=FALSE]	# missing rate control at 0.05.
		 # ---2, monomorphic site control---#
		 singleSNV_MAF <- iter[,colMeans(.SD,na.rm=TRUE),.SDcols = -'rn']
		 # iter <- iter[,.SD,.SDcols = setdiff(colnames(iter), 'rn')]
		 iter <- iter[,!names(singleSNV_MAF)[singleSNV_MAF == 0], with=FALSE]
		 # ---3, we skip those genes with 0 or 1 SNVs---# 
		 if ( ncol(iter) <= 2) return()	# we skip those genes with 0 or 1 SNVs. because we have 'patientID', so the cutoff is set at <= 2.
		 
		 # ---setkey for iter obj---#
		 setkeyv(iter,"rn")
		 # #----------------------#
		 # data <- data.frame(id = 1:nrow(Y),Y = Y, X = X,check.names=F)
		 # data_long <- plyr::arrange(reshape(data,idvar = "id",varying = grep('Y\\.',colnames(data),perl=T),direction = "long",v.names="Y",times = 1:longi_n), id )
		 data_long <- iter[phe_cov_longitudinal_longFormat]	# A JOIN here.
		 # stopifnot(  colnames(data_long)[grep("^exm",colnames(data_long))] %in% colnames(iter[,!"rn",with=FALSE]))	# it's not robust, if the SNV is not start with '^exm', so only check for a few debugging samples to test the algorithm is right.
		 # we already have data_long ready for real data.
		 #-----if want to use geescore_for-----#
		 data_long <- mutate(data_long, intercept=1)	# this is the cov.
		 # data_long <- tbl_dt(data_long)
		 stopifnot( covNames %in% colnames(data_long) )
		 parNames <- colnames(iter[,!"rn",with=FALSE])
		 #
		 
		 
		 #---run geescore_for---------------------------------------#
		 system.time( geescore_Result <- geescore_for(nCluster = NULL,data = data_long, id = "rn", y = responseVarName, 
					cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian", .parallel=parallelFlag, allowCovNA = T, allowParNA = T) )	# change the y here if you have multiple longitudinal responses.

			#---------------------------#
			# after we got the geescore_Result from "Test_simulation_data_execute_RV_version.r"
			pSSU <- SumSqU(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pSSUw <- SumSqUw(U = geescore_Result$U, CovS = geescore_Result$Cov)
			# pScore <- Score(U = geescore_Result$U, CovS = geescore_Result$Cov)
			#---modified pScore.
			pScore <- tryCatch( Score(U = geescore_Result$U, CovS = geescore_Result$Cov),
								# warning = identity, 
								error = function(x) NaN	# omit error msg, just output NaN. So it can still be rbinded with other numeric pvalues.						
								)	# this score version do some normalization on too small values in Cov matrix.
			pSum <- Sum(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pUminP <- UminP(U = geescore_Result$U, CovS = geescore_Result$Cov)
		 
			return( c(pSSU =pSSU, pSSUw =pSSUw, pScore=pScore, pSum=pSum, pUminP=pUminP,  geescoreP = geescore_Result$geeScoreTest_out$pval) )
		 
		}
		names(testResult) <- allGeneNames	# stored by gene name for each list element.
	

} else {	# include aSPU family tests
# #---depending on parallel_scheme, setup the parallelFlag-----------------#
# # if we have multiple nodes available, we can use another parallel layer for 'geescore_for' computing (within one slave node for one gene); if not, use serial here.
# # similar situation for SPU3_withScoGee_RV_version.
# if (exists("MultipleNodes_profile")) {	# we have multipleNodes ready.
	# parallelFlag <- TRUE 
# } else {
	# parallelFlag <- FALSE
# }
	if (!parallel_over_gene) {	
		parallelFlag <- TRUE
		#--run---#
		# testResult <- foreach(iter = iterObj,.errorhandling = "pass",.verbose=T) %dopar% {
		testResult <- foreach(iter = iterObj,.errorhandling = "pass",.verbose=T) %do% {	# gene level in serial, within each gene, the computation is paralleled.

		 # registerDoMC(detectCores())	# for each in cl_mother
		 # iter = nextElem(iterObj)	# for debugging usage.
		 #----some limited QC here, more extensive QC should be done outside of this program----#
		 # ---1, NA control----#
		 NaRatios <- apply(iter,2,function(x) mean(is.na(x)) )
		 iter <- iter[,NaRatios <= 0.05,with=FALSE]	# missing rate control at 0.05.
		 # ---2, monomorphic site control---#
		 singleSNV_MAF <- iter[,colMeans(.SD,na.rm=TRUE),.SDcols = -'rn']
		 # iter <- iter[,.SD,.SDcols = setdiff(colnames(iter), 'rn')]
		 iter <- iter[,!names(singleSNV_MAF)[singleSNV_MAF == 0], with=FALSE]
		 # ---3, we skip those genes with 0 or 1 SNVs---# 
		 if ( ncol(iter) <= 2) return()	# we skip those genes with 0 or 1 SNVs. because we have 'patientID', so the cutoff is set at <= 2.
		 
		 # ---setkey for iter obj---#
		 setkeyv(iter,"rn")
		 # #----------------------#
		 # data <- data.frame(id = 1:nrow(Y),Y = Y, X = X,check.names=F)
		 # data_long <- plyr::arrange(reshape(data,idvar = "id",varying = grep('Y\\.',colnames(data),perl=T),direction = "long",v.names="Y",times = 1:longi_n), id )
		 data_long <- iter[phe_cov_longitudinal_longFormat]	# A JOIN here.
		 # stopifnot(  colnames(data_long)[grep("^exm",colnames(data_long))] %in% colnames(iter[,!"rn",with=FALSE]))	# it's not robust, if the SNV is not start with '^exm', so only check for a few debugging samples to test the algorithm is right.
		 # we already have data_long ready for real data.
		 #-----if want to use geescore_for-----#
		 data_long <- mutate(data_long, intercept=1)	# this is the cov.
		 # data_long <- tbl_dt(data_long)
		 stopifnot( covNames %in% colnames(data_long) )
		 parNames <- colnames(iter[,!"rn",with=FALSE])
		 #
		 
		 
		 #---run geescore_for---------------------------------------#
		 system.time( geescore_Result <- geescore_for(nCluster = NULL,data = data_long, id = "rn", y = responseVarName, 
					cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian", .parallel=parallelFlag, allowCovNA = T, allowParNA = T) )	# change the y here if you have multiple longitudinal responses.

			#---------------------------#
			# after we got the geescore_Result from "Test_simulation_data_execute_RV_version.r"
			pSSU <- SumSqU(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pSSUw <- SumSqUw(U = geescore_Result$U, CovS = geescore_Result$Cov)
			# pScore <- Score(U = geescore_Result$U, CovS = geescore_Result$Cov)
			#---modified pScore.
			pScore <- tryCatch( Score(U = geescore_Result$U, CovS = geescore_Result$Cov),
								# warning = identity, 
								error = function(x) NaN	# omit error msg, just output NaN. So it can still be rbinded with other numeric pvalues.						
								)	# this score version do some normalization on too small values in Cov matrix.
			pSum <- Sum(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pUminP <- UminP(U = geescore_Result$U, CovS = geescore_Result$Cov)

			 #----------if usePermutedU, we will generate null U empirical distribution using permutation strategy--------#
			if (usePermutedU) {
			 system.time(
			 null_Ustatistics <- RV_permute_test(nCluster = NULL,data = data_long, id = "rn", y = responseVarName, 
						# cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian",N_Perm=b_times, .parallel=FALSE, .seed=1234, cl = cl)	# this will be very slow without parallel.
						cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian",N_Perm=b_times, .parallel=parallelFlag, .seed=seed, cl = cl)	# need parallel here.
			 )
			}

			#-----------aSPU test here-----------------#
			if (usePermutedU) {
			 system.time(spu_Result<-SPU3_withScoGee_RV_version(U=geescore_Result$U, V=geescore_Result$Cov, gamma=c(0:8,Inf),B=b_times, .seed=seed, usePermutationNullU = TRUE, permutationNullU = null_Ustatistics, .parallel_overGamma = parallelFlag, .parallel_overB = parallelFlag) )	# for parallel within gene scheme, we will parallel over gamma and B to fully utilize the cpus resources
			} else {
			 system.time(spu_Result<-SPU3_withScoGee_RV_version(U=geescore_Result$U, V=geescore_Result$Cov, gamma=c(0:8,Inf),B=b_times, .seed=seed, usePermutationNullU = FALSE, .parallel_overGamma = parallelFlag, .parallel_overB = parallelFlag) )
			}

		 
		# return( c(pSSU =pSSU,pSSUw =pSSUw,pScore=pScore,pSum =pSum,pUminP=pUminP
			# ,geescoreP = geescore_Result$geeScoreTest_out$pval,
			# UminP = spu_Result$minP, aveP= spu_Result$spu[['1']],  aveP_weighted = spu_Result$spu.w[['1']], KM_P = spu_Result$spu[['2']], geescoreP_withDiagSigma = spu_Result$spu.w[['2']], spu_Inf_P = spu_Result$spu[['Inf']], spu_weighted_Inf_P =  spu_Result$spu.w[['Inf']], aspu_P = spu_Result$spu[['aspu']], aspu_weighted_P = spu_Result$spu.w[['aspu.w']], aspu_sco_P = spu_Result$spu.sco[['aspu.sco']],aspu_weighted_sco_P = spu_Result$spu.w.sco[['aspu.w.sco']], aspu_aspuw.sco_P = spu_Result$Paspu_aspuw.sco, aspu_aspuw_P = spu_Result$Paspu_aspuw   ) )
		return( c(pSSU =pSSU, pSSUw =pSSUw, pScore=pScore, pSum=pSum, pUminP=pUminP
			, geescoreP = geescore_Result$geeScoreTest_out$pval, UminP = spu_Result$minP, spu1= spu_Result$spu[['1']],  spu1w = spu_Result$spu.w[['1']], spu2 = spu_Result$spu[['2']], spu2w = spu_Result$spu.w[['2']], spu_Inf = spu_Result$spu[['Inf']], spu_Inf_w =  spu_Result$spu.w[['Inf']], aspu_P = spu_Result$spu[['aspu']], aspu_w_P = spu_Result$spu.w[['aspu.w']], aspu_sco_P = spu_Result$spu.sco[['aspu.sco']],aspu_w_sco_P = spu_Result$spu.w.sco[['aspu.w.sco']], aspu_aspuw.sco_P = spu_Result$Paspu_aspuw.sco, aspu_aspuw_P = spu_Result$Paspu_aspuw   ) )
		 
		}
		names(testResult) <- allGeneNames	# stored by gene name for each list element.
		
	} else {
		parallelFlag <- FALSE
			#--run---#
		testResult <- foreach(iter = iterObj,.errorhandling = "pass",.verbose=T) %dopar% {
		# testResult <- foreach(iter = iterObj,.errorhandling = "pass",.verbose=T) %do% {	# gene level in serial, within each gene, the computation is paralleled.

		 # registerDoMC(detectCores())	# for each in cl_mother
		 # iter = nextElem(iterObj)	# for debugging usage.
		 #----some limited QC here, more extensive QC should be done outside of this program----#
		 # ---1, NA control----#
		 NaRatios <- apply(iter,2,function(x) mean(is.na(x)) )
		 iter <- iter[,NaRatios <= 0.05,with=FALSE]	# missing rate control at 0.05.
		 # ---2, monomorphic site control---#
		 singleSNV_MAF <- iter[,colMeans(.SD,na.rm=TRUE),.SDcols = -'rn']
		 # iter <- iter[,.SD,.SDcols = setdiff(colnames(iter), 'rn')]
		 iter <- iter[,!names(singleSNV_MAF)[singleSNV_MAF == 0], with=FALSE]
		 # ---3, we skip those genes with 0 or 1 SNVs---# 
		 if ( ncol(iter) <= 2) return()	# we skip those genes with 0 or 1 SNVs. because we have 'patientID', so the cutoff is set at <= 2.
		 
		 # ---setkey for iter obj---#
		 setkeyv(iter,"rn")
		 # #----------------------#
		 # data <- data.frame(id = 1:nrow(Y),Y = Y, X = X,check.names=F)
		 # data_long <- plyr::arrange(reshape(data,idvar = "id",varying = grep('Y\\.',colnames(data),perl=T),direction = "long",v.names="Y",times = 1:longi_n), id )
		 data_long <- iter[phe_cov_longitudinal_longFormat]	# A JOIN here.
		 # stopifnot(  colnames(data_long)[grep("^exm",colnames(data_long))] %in% colnames(iter[,!"rn",with=FALSE]))	# it's not robust, if the SNV is not start with '^exm', so only check for a few debugging samples to test the algorithm is right.
		 # we already have data_long ready for real data.
		 #-----if want to use geescore_for-----#
		 data_long <- mutate(data_long, intercept=1)	# this is the cov.
		 # data_long <- tbl_dt(data_long)
		 stopifnot( covNames %in% colnames(data_long) )
		 parNames <- colnames(iter[,!"rn",with=FALSE])
		 #
		 
		 
		 #---run geescore_for---------------------------------------#
		 system.time( geescore_Result <- geescore_for(nCluster = NULL,data = data_long, id = "rn", y = responseVarName, 
					cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian", .parallel=parallelFlag, allowCovNA = T, allowParNA = T) )	# change the y here if you have multiple longitudinal responses.

			#---------------------------#
			# after we got the geescore_Result from "Test_simulation_data_execute_RV_version.r"
			pSSU <- SumSqU(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pSSUw <- SumSqUw(U = geescore_Result$U, CovS = geescore_Result$Cov)
			# pScore <- Score(U = geescore_Result$U, CovS = geescore_Result$Cov)
			#---modified pScore.
			pScore <- tryCatch( Score(U = geescore_Result$U, CovS = geescore_Result$Cov),
								# warning = identity, 
								error = function(x) NaN	# omit error msg, just output NaN. So it can still be rbinded with other numeric pvalues.						
								)	# this score version do some normalization on too small values in Cov matrix.
			pSum <- Sum(U = geescore_Result$U, CovS = geescore_Result$Cov)
			pUminP <- UminP(U = geescore_Result$U, CovS = geescore_Result$Cov)

			 #----------if usePermutedU, we will generate null U empirical distribution using permutation strategy--------#
			if (usePermutedU) {
			 system.time(
			 null_Ustatistics <- RV_permute_test(nCluster = NULL,data = data_long, id = "rn", y = responseVarName, 
						# cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian",N_Perm=b_times, .parallel=FALSE, .seed=1234, cl = cl)	# this will be very slow without parallel.
						cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian",N_Perm=b_times, .parallel=parallelFlag, .seed=seed, cl = cl)	# need parallel here.
			 )
			}

			#-----------aSPU test here-----------------#
			if (usePermutedU) {
			 system.time(spu_Result<-SPU3_withScoGee_RV_version(U=geescore_Result$U, V=geescore_Result$Cov, gamma=c(0:8,Inf),B=b_times, .seed=seed, usePermutationNullU = TRUE, permutationNullU = null_Ustatistics, .parallel_overGamma = parallelFlag, .parallel_overB = FALSE) )
			} else {
			 system.time(spu_Result<-SPU3_withScoGee_RV_version(U=geescore_Result$U, V=geescore_Result$Cov, gamma=c(0:8,Inf),B=b_times, .seed=seed, usePermutationNullU = FALSE, .parallel_overGamma = parallelFlag, .parallel_overB = FALSE) )
			}

		 
		# return( c(pSSU =pSSU,pSSUw =pSSUw,pScore=pScore,pSum =pSum,pUminP=pUminP
			# ,geescoreP = geescore_Result$geeScoreTest_out$pval,
			# UminP = spu_Result$minP, aveP= spu_Result$spu[['1']],  aveP_weighted = spu_Result$spu.w[['1']], KM_P = spu_Result$spu[['2']], geescoreP_withDiagSigma = spu_Result$spu.w[['2']], spu_Inf_P = spu_Result$spu[['Inf']], spu_weighted_Inf_P =  spu_Result$spu.w[['Inf']], aspu_P = spu_Result$spu[['aspu']], aspu_weighted_P = spu_Result$spu.w[['aspu.w']], aspu_sco_P = spu_Result$spu.sco[['aspu.sco']],aspu_weighted_sco_P = spu_Result$spu.w.sco[['aspu.w.sco']], aspu_aspuw.sco_P = spu_Result$Paspu_aspuw.sco, aspu_aspuw_P = spu_Result$Paspu_aspuw   ) )
		return( c(pSSU =pSSU, pSSUw =pSSUw, pScore=pScore, pSum=pSum, pUminP=pUminP
			, geescoreP = geescore_Result$geeScoreTest_out$pval, UminP = spu_Result$minP, spu1= spu_Result$spu[['1']],  spu1w = spu_Result$spu.w[['1']], spu2 = spu_Result$spu[['2']], spu2w = spu_Result$spu.w[['2']], spu_Inf = spu_Result$spu[['Inf']], spu_Inf_w =  spu_Result$spu.w[['Inf']], aspu_P = spu_Result$spu[['aspu']], aspu_w_P = spu_Result$spu.w[['aspu.w']], aspu_sco_P = spu_Result$spu.sco[['aspu.sco']],aspu_w_sco_P = spu_Result$spu.w.sco[['aspu.w.sco']], aspu_aspuw.sco_P = spu_Result$Paspu_aspuw.sco, aspu_aspuw_P = spu_Result$Paspu_aspuw   ) )
		 
		}
		names(testResult) <- allGeneNames	# stored by gene name for each list element.

	}



}






#----save the result---#
setwd( sprintf('%s/%s', programRootPWD, 'OutputData') )
filename_essence <- sub("\\.Rdata$","",filename_genotype,perl=T,ignore.case = TRUE)
# 
if (usePermutedU) {
	Which_U_string <- "usePermutationU"
} else {
	Which_U_string <- "useSimulationU"
}

save(testResult, file = sprintf("%s_response:%s_%s_wkCor:%s_bTime:%d_%s.Rdata",prefix,responseVarName,filename_essence,workingCorstr_input,b_times,Which_U_string ))

#---safely end up the virtual cluster---#
if (exists("cl")) stopCluster(cl)
cat("Finished at ", date(),"\n" )
cat("=====================================================\n")

# ##################################################################################
				  