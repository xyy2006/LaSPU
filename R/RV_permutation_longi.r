## Content: do permutation to get the uNull.
##			1st, complement the data if has missings;
##			2nd, permute y across individuals while holding the within-subject measurements in the same order (e.g. timepoints 1~4.)
##			3rd, use 'geescore_for' to test permuted dataset.
###############################################################################
# usage: source this script in main script.
###############################################################################
RV_permute_test <- function(nCluster = NULL,data = data_long, id = "id", y = "Y", 
			cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian",N_Perm=1000, .parallel=F, cl = NULL,.seed=1234) {
	# browser()
	# > match.call(RV_permute_test,sys.call())
	# RV_permute_test(nCluster = NULL, data = data_long, id = "id",
	# y = "Y", cov = covNames, par = parNames, waves = "time",
	# corstr = workingCorstr_input, biasCorrection = T, commonVariance = T,
	# family = "gaussian", N_Perm = b_times, .parallel = T, .seed = 1234)
	# data <- data[sample(nrow(data),size = 0.8 * nrow(data) )]	# simulate missing data, for debugging.
	#-----------------------#
	# old obsolete
	# if (id == "id") {
		# data <- plyr::arrange(data, id, eval(as.name(waves)) )	# plyr allow expr in ...
	# } else 	data <- plyr::arrange(data, eval(as.name(id)), eval(as.name(waves)) )
	#-----------------------#
	if (!.parallel) `%dopar%` <- `%do%`
	stopifnot(is.data.frame(data))	# we need data to be data.frame instead of matrix.
	data <- tbl_dt(data)
	arrange_call <- (sprintf("plyr::arrange(data, %s, %s)", id, waves))	# arrange by id and then waves.
	data <- eval(parse(text = arrange_call))
	#
	if (is.null(nCluster)) nCluster <- length(unique(data[[eval(id)]]))	# for data.frame data, how many nCluster we have after all data filtering.
	#
	waves_by_id <- table(data[[eval(id)]])
	max_longiN <- max(waves_by_id)
	#
	if (all (waves_by_id == max(waves_by_id) )) {
		full_obs_flag <- TRUE	# perfect, we have full obs data, which will facilitate the below computation.
	} else full_obs_flag <- FALSE	

	#
	if (!full_obs_flag){	# only do it when we have missing data.
		collect_all_measurements_seqs <- function(data_chunk = subset(data,subset = eval(as.name(id)) == 1) ){
		  # browser()
		  measurements_seq <- data_chunk[[waves]]
		  return(measurements_seq)
		}
		if (id == "id") {
			measurements_seq_pool <- dlply(.data = data, .variables =.(id), .fun =collect_all_measurements_seqs, .progress = "text", .parallel = FALSE, .paropts = list(.packages=c("MASS","gdata","plyr")), .inform=F)
		} else {
			measurements_seq_pool <- dlply(.data = data, .variables =.(eval(as.name(id))), .fun =collect_all_measurements_seqs, .progress = "text", .parallel = FALSE, .paropts = list(.packages=c("MASS","gdata","plyr")),.inform=F)
		}
	

		# measurements_seq_pool	# will be used later in building subject specific Rw.
		measurements_seq_pool_pasted <- sapply(measurements_seq_pool,paste,collapse='')	# for N subjects, length = n
		measurements_seq_categories <- names( table(measurements_seq_pool_pasted) )
		#
		
		#---define the max timepoints we can get:---#
		max_groupName <- names(which.max(sapply(measurements_seq_categories,nchar)))	# "1234"
		max_groupName_vec <- strsplit(max_groupName,split="")[[1]]
		#---get the subjectID in groups by measurements_seq_categories.---#
		subjectID_group <- foreach ( measurement_seq = measurements_seq_categories) %dopar% {
			names(measurements_seq_pool_pasted)[ which(measurements_seq_pool_pasted %in% measurement_seq) ]
		}
		names(subjectID_group) <- measurements_seq_categories
		
		#---complement missing data to make 'data' complete dimension---#
		# then permute can be easier for 'full' obs.
		#---a subfunc to complement missing data for each obs who has missings.---#
		complement_missingObs <- function(ids_in_a_group = subjectID_group[[2]], groupNames = names(subjectID_group)[2] ){
			# browser()
			if (groupNames == max_groupName) return();	# complete obs dont need complement any more.
			#------------------#
			N_timepoints <- nchar(groupNames)
			groupNames_vec <- strsplit(groupNames,split="")[[1]]
			filter_call <- sprintf("dplyr::filter(data, %s %s ids_in_a_group)",id, "%in%")
			data_ids_in_a_group <- eval(parse(text=filter_call))	# id and wave are in good order due to the arrange function used in the beginning.
			#
			stopifnot( (data_ids_in_a_group[,id,with=F]) == rep(ids_in_a_group,each = N_timepoints) )	# check if the order is neat and the same.
			#
			#---start complement---#
			N_timepoints_complement <- max_longiN - N_timepoints
			groupNames_complement_vec <- setdiff(max_groupName_vec,groupNames_vec)	# the complement time points. For example, if we have only 1 and 2, then 3 and 4 are complement time points to 1:4.
			groupNames_complement_vec_int <- as.integer(groupNames_complement_vec)
			unique_id_row <- !duplicated(data_ids_in_a_group[,id,with=F])
			data_ids_in_a_group_unique_row <- data_ids_in_a_group[unique_id_row]
			#
			data_ids_in_a_group_complemented <- c()
			for (i in seq(N_timepoints_complement) ) {
				# data_ids_in_a_group_unique_row[,`:=`(waves= groupNames_complement_vec_int[i], y=NA),with=F]	# cannot work with string name inside a name symbol.
				data_ids_in_a_group_unique_row[,(waves):= groupNames_complement_vec_int[i],with=F]	# give the complement time point value
				data_ids_in_a_group_unique_row[,(y):= NA,with=F]	# since it's missing value, we should give y a value of NA to keep format nice.
				data_ids_in_a_group_complemented <- rbind(data_ids_in_a_group_complemented,data_ids_in_a_group_unique_row)
			} 
			return(data_ids_in_a_group_complemented)	# return only the complemented one, not the rbinded complete one.
		}
	
		#
		# Browse[1]> names(subjectID_group)
		# [1] "1"    "12"   "123"  "1234" "124"  "13"   "134"  "14"   "2"    "23"
		# [11] "234"  "24"   "3"    "34"   "4"
		#---execute the complement process---#
		data_complemented <- mapply(complement_missingObs,subjectID_group,names(subjectID_group) )	# the complemented data list with names
		#
		data_expanded_to_fullDim <- rbind(data, do.call(rbind, data_complemented) )
		stopifnot( all(data_expanded_to_fullDim$id %>>% table == max_longiN) )
		data <- data_expanded_to_fullDim	# update the data obj.

	}
	arrange_call <- (sprintf("plyr::arrange(data, %s, %s)", id, waves))	# arrange by id and then waves.
	data <- eval(parse(text = arrange_call))	# re-arrange to make sure everything is neat and clean.
	subjectID <- unlist( data[,id,with=F], use.names=FALSE) %>>% unique
	# print(subjectID)
	data_backup <- copy(data)	# save this one, so data can be changed in any way.
	# #---a subfunc to do permutation---#
	# permute_within_group <- function(ids = data[,id,with=F] ) {
		# sample(ids)	# return the permuted new group of ids.
	# }
	# return(list(data_expanded = data, subjectID_expanded = subjectID))
	# browser()
	# #-------start execute permutation to generate permutedNull_Ustatistics---#
	# registerDoRNG(.seed)	# guarantee parallel reproducible permutations.
	# clusterExport(cl,list = ls(),envir = environment())
	# clusterExport(cl,list = ls())
	if (.parallel & !is.null(cl) ) {
	  clusterEvalQ(cl, source("head.r"))	# may change in the future.
	} else if (.parallel & is.null(cl) ) {
	  registerDoMC(detectCores())
	}
	permutedNull_Ustatistics <- foreach(i = seq(N_Perm), .combine = rbind,.verbose=T,.export = "geescore_for",.inorder=FALSE) %dopar% {	# order doesn't matter. all other objects in the calling env will be exported except 'geescore_for', which is in the parent env of calling env, so need to import by .export arg.
		subjectID_permuted <- sample( subjectID )
		#
		#---we can start mapping between original and permuted---#
		ids_permuted_vec <- paste(rep(subjectID_permuted,each = max_longiN),seq(max_longiN),sep='_')
		ids_vec <- paste(rep(subjectID,each = max_longiN),seq(max_longiN),sep='_')

		index <- match( ids_permuted_vec, ids_vec )
		#-----------permute y+cov or par-------------------------#
		# Key: we want to hold the y+covariate part permuted together, so either permute them across subjects, or simply permute the par(SNP,SNV) across subjects. Keep the longitudinal order.
		#
		# new_Y_vec <- data_backup[index,y,with=F]	# the new y in 
		# data[,y:=new_Y_vec,with=F]	# just change y part, keep other parts intact is our strategy.
		 # # object.size(data)
		# # [1] 4.6 MB
		new_par_vec <- data_backup[index,par,with=F]	# the new par in 
		data[,par:=new_par_vec,with=F]	# just change par part, keep other parts intact is our strategy. so it doesn't break the relationship between y and covariates if any.
		 # object.size(data)
		# [1] 4.6 MB
		
		#-----------start do geescore_for on permuted data-------#
		geescore_Result <- geescore_for(nCluster = NULL,data = data, id = id, y = y, cov = cov , par = par,waves = waves, corstr=corstr, biasCorrection=biasCorrection, commonVariance=commonVariance, family=family, .parallel=FALSE)	# set .parallel to false is very important. 
		 # obj_size(geescore_Result)
		# [1] 26 kB
		return( t(geescore_Result$U) )
		# t(geescore_Result$U)
		# X.1       X.2      X.3      X.4      X.5       X.6       X.7
		# [1,] -18.40682 -11.98736 23.02139 3.387642 -16.2395 0.8987657 -2.580056
		# X.8       X.9      X.10      X.11      X.12     X.13     X.14
		# [1,] 17.06902 -21.53267 -55.56688 -35.59789 -1.438798 11.85432 29.67662

	}
	
	return(permutedNull_Ustatistics)
}
	
 



#==============================================================================#
# original version, permutate within group. (e.g. timepoints: 12, 1234, 234, 1)
# RV_permute <- function(nCluster = NULL,data = data_long, id = "id", y = "Y", cov = covNames , par = parNames,waves = "time", corstr=workingCorstr_input, biasCorrection=T, commonVariance=T, family="gaussian", N_Perm = 1000, .parallel= F) {
	# # browser()
	# data <- data[sample(nrow(data),size = 0.8 * nrow(data) )]	# simulate missing data
	# #-----------------------#
	# # old obsolete
	# # if (id == "id") {
		# # data <- plyr::arrange(data, id, eval(as.name(waves)) )	# plyr allow expr in ...
	# # } else 	data <- plyr::arrange(data, eval(as.name(id)), eval(as.name(waves)) )
	# #-----------------------#
	# data <- tbl_dt(data)
	# arrange_call <- (sprintf("plyr::arrange(data, %s, %s)", id, waves))	# arrange by id and then waves.
	# data <- eval(parse(text = arrange_call))
	# #
	# if (is.null(nCluster)) nCluster <- length(unique(data[[eval(id)]]))	# for data.frame data, how many nCluster we have after all data filtering.
	# #
	# waves_by_id <- table(data[[eval(id)]])
	# max_longiN <- max(waves_by_id)
	# #
	# if (all (waves_by_id == max(waves_by_id) )) {
		# full_obs_flag <- TRUE	# perfect, we have full obs data, which will facilitate the below computation.
	# } else full_obs_flag <- FALSE	

	# #
	# if (!full_obs_flag){	# only do it when we have missing data.
		# collect_all_measurements_seqs <- function(data_chunk = subset(data,subset = eval(as.name(id)) == 1) ){
		  # # browser()
		  # measurements_seq <- data_chunk[[waves]]
		  # return(measurements_seq)
		# }
		# if (id == "id") {
			# measurements_seq_pool <- dlply(.data = data, .variables =.(id), .fun =collect_all_measurements_seqs, .progress = "text", .parallel = FALSE, .paropts = list(.packages=c("MASS","gdata","plyr")), .inform=F)
		# } else {
			# measurements_seq_pool <- dlply(.data = data, .variables =.(eval(as.name(id))), .fun =collect_all_measurements_seqs, .progress = "text", .parallel = FALSE, .paropts = list(.packages=c("MASS","gdata","plyr")),.inform=F)
		# }
	# }

	# # measurements_seq_pool	# will be used later in building subject specific Rw.
	# measurements_seq_pool_pasted <- sapply(measurements_seq_pool,paste,collapse='')	# for N subjects, length = n
	# measurements_seq_categories <- names( table(measurements_seq_pool_pasted) )
	
	# #---get the subjectID in groups by measurements_seq_categories.---#
	# subjectID_group <- foreach ( measurement_seq = measurements_seq_categories) %dopar% {
		# names(measurements_seq_pool_pasted)[ which(measurements_seq_pool_pasted %in% measurement_seq) ]
	# }
	# names(subjectID_group) <- measurements_seq_categories
	
	# browser()
	# #---a subfunc to do permutation within each group---#
	# permute_within_group <- function(groupList = subjectID_group[[1]]) {
		# # browser()
		# sample(groupList)	# return the permuted new group of ids.
		# # Browse[3]> sample(groupList)
		# # [1] 2332  212  464  860 2404 1557 2546 1875 1479 1674 1556 1821 2329 2747  932
		# # [16] 2146 2629 1601  122 2362 2935	
	# }
	
	# #---start to permute---#
	# result <- times(N_Perm) %dopar% {
		# subjectID_group_permuted <- sapply(subjectID_group, permute_within_group)
		# stopifnot( sapply(subjectID_group_permuted, length) == sapply(subjectID_group, length) )	# after permute, each group should have same length.
		# #
		# foreach ( ids_in_a_group = subjectID_group, ids_in_permuted_group = subjectID_group_permuted, groupNames =  names(subjectID_group) ) %do% {
			# N_timepoints <- nchar(groupNames)
			# filter_call <- sprintf("dplyr::filter(data, %s %s ids_in_a_group)",id, "%in%")
			# data_ids_in_a_group <- eval(parse(text=filter_call))	# id and wave are in good order due to the arrange function used in the beginning.
			# #
			# stopifnot(data_ids_in_a_group[,id,with=F] == rep(ids_in_a_group,each = N_timepoints) )	# check if the order is neat and the same.
			# stopifnot(ids_in_a_group %in% ids_in_permuted_group )
			# #---we can start mapping between original and permuted---#
			# ids_in_permuted_group_vec <- paste(rep(ids_in_permuted_group,each = N_timepoints),seq(N_timepoints),sep='_')
			# ids_in_a_group_vec <- paste(rep(ids_in_a_group,each = N_timepoints),seq(N_timepoints),sep='_')
			
			# index <- match( ids_in_permuted_group_vec, ids_in_a_group_vec )
			# #
			# stopifnot(data_ids_in_a_group[index,id,with=F] == rep(ids_in_permuted_group,each = N_timepoints))	# check the order in permuted is neat and the same.
			# new_Y_vec <- data_ids_in_a_group[index,y,with=F]	# the new y in 
			# data_ids_in_a_group[,y:=new_Y_vec,with=F]	# just change y part, keep other parts intact is our strategy.
			# # Browse[1]> data_ids_in_a_group[,y,with=F]
					 # # Y
			# # 1: 0.3149450
			# # 2: 0.8106353
			# # 3: 4.0606624
			# # 4: 2.3290820
			# # 5: 2.3407171
			# # ---
			# # 154: 2.1826904
			# # 155: 2.3076043
			# # 156: 2.5810576
			# # 157: 0.6251735
			# # 158: 0.4411334
			# #####
			# # change to
			# #####
			# # Browse[1]> data_ids_in_a_group[,y,with=F]
			# # Y
			# # 1: 3.712802
			# # 2: 2.880446
			# # 3: 1.746102
			# # 4: 1.388098
			# # 5: 1.778007
			# # ---
			# # 154: 3.763420
			# # 155: 2.147173
			# # 156: 2.614559
			# # 157: 3.857099
			# # 158: 3.573140


		# }
		
	# }
	
 

# }
