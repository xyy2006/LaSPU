#!/bin/bash
# the Bash shell script wrapper for `ModelComputation_main` program.
#---find the log file we just generated---#
# echo $1 $2 ${10}  # test print several args.
./assoc-aSPU.r $@ |& tee temp.log
logfile=(`ls -ahlrt *.log|grep $USER | tail -n 2 |awk '{print $9}'`)
for index in ${!logfile[@]}; 
do
  value=${logfile[$index]}
  if [[ "${value}" != "temp.log" ]]
  then
    index_return=${index}
  fi
done
cat temp.log > ${logfile[$index_return]} 
unlink temp.log
echo "========================================================"
echo "Log file name saved to: ./"${logfile[$index_return]}
# echo "Print it on screen now:"
# cat ${logfile[$index_return]} 

# ./ModelComputation_main $@ >(tee stdout.log) 2> >(tee stderr.log >&2) # split the output and erroutput.
# #---deal with the logfile just generated by the pipeline---# since it is no longer necessary with `tee` cmd above.
# logfile=`ls -ahlrt *.log|grep $USER | tail -n 1 |awk '{print $9}'`
# echo "Log file name is:" $logfile
# echo "Print it on screen now:"
# cat ${logfile}