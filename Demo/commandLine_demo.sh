##################################################################################################################
# Content: introduce how to use the program, e.g. commandline operation demo
# Author: Yang Yang
##################################################################################################################

# ./OutputData/singleNode_simulation_response:trigs_longi_genotpe_data_sample_wkCor:independence_bTime:1000_useSimulationU.Rdata
#---1, single node parallel run using simulation method based empirical null U distribution---#
./assoc-aSPU.r --filename_genotype genotpe_data_sample.Rdata \
               --filename_annotation genotype_annotation.Rdata \
               --filename_phenotype phe_cov_longitudinal_longFormat_sampleData.Rdata \
               --responseVarName trigs_longi \
               --covariateVarFile covariateNames.txt \
               --workingCorstr_input independence \
               --include_aSPU TRUE \
               --b_times 1e3 \
               --seed 123456 \
               --usePermutedU FALSE \
               --parallel_scheme SNOW \
               --parallel_over_gene TRUE \
               --jobExecuteOnMultipleNodes FALSE \
               --prefix singleNode_simulation

							 
							 
# ./OutputData/singleNode_without_aSPU_response:trigs_longi_genotpe_data_sample_wkCor:independence_bTime:1000_useSimulationU.Rdata
#---2, single node parallel run with 'include_aSPU' = FALSE ---#
# only report default 5 tests results
./assoc-aSPU.r --filename_genotype genotpe_data_sample.Rdata \
               --filename_annotation genotype_annotation.Rdata \
               --filename_phenotype phe_cov_longitudinal_longFormat_sampleData.Rdata \
               --responseVarName trigs_longi \
               --covariateVarFile covariateNames.txt \
               --workingCorstr_input independence \
               --include_aSPU FALSE \
               --b_times 1e3 \
               --seed 123456 \
               --usePermutedU FALSE \
               --parallel_scheme SNOW \
               --parallel_over_gene TRUE \
               --jobExecuteOnMultipleNodes FALSE \
               --prefix singleNode_without_aSPU
							 
# no result							 
#---3, single nodes parallel run using permutation method based empirical null U distribution---#
./assoc-aSPU.r --filename_genotype genotpe_data_sample.Rdata \
               --filename_annotation genotype_annotation.Rdata \
               --filename_phenotype phe_cov_longitudinal_longFormat_sampleData.Rdata \
               --responseVarName trigs_longi \
               --covariateVarFile covariateNames.txt \
               --workingCorstr_input independence \
               --include_aSPU TRUE \
               --b_times 1e3 \
               --seed 123456 \
               --usePermutedU TRUE \
               --parallel_scheme SNOW \
               --parallel_over_gene FALSE \
               --jobExecuteOnMultipleNodes FALSE \
               --prefix singleNode_permutation
# no result
#---4, multiple nodes parallel run using permutation method based empirical null U distribution---#
./assoc-aSPU.r --filename_genotype genotpe_data_sample.Rdata \
               --filename_annotation genotype_annotation.Rdata \
               --filename_phenotype phe_cov_longitudinal_longFormat_sampleData.Rdata \
               --responseVarName trigs_longi \
               --covariateVarFile covariateNames.txt \
               --workingCorstr_input independence \
               --include_aSPU TRUE \
               --b_times 1e3 \
               --seed 123456 \
               --usePermutedU TRUE \
               --parallel_scheme SNOW \
               --parallel_over_gene FALSE \
               --jobExecuteOnMultipleNodes TRUE \
               --MultipleNodes_profile MultipleNodes_profile2.txt \
               --prefix multipleNodes_permutation
							 