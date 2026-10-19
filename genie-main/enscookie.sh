#!/bin/bash
#
# *** SUBMIT ENSEMBLE ***
#
# NOTE: order of ensemble member extensions is: 
#       o n m

# specify base-config filename
BASECONFIG=''
# specify path to any user-configs ('/' otherwise)
USERCONFIGPATH=''
# specify experiment (user-config) filename, excluding the .xy ensemble member extension
ENSEMBLEID=''
# specify run duration (integer years)
YEARS=''
# specify any restart name (empty string otherwise)
RESTARTID=''
# set x-axis parameter	== max value in decimal (1-15) of THIRD hex digit
mmax=1
# set y-axis parameter 	== max value in decimal (1-15) of SECOND hex digit
nmax=1
# set z-axis parameter	== max value in decimal (1-15) of FIRST hex digit
omax=1

# initialize loop counter m
o=1
while [ $o -le $omax ]; do
  # initialize loop counter n
  n=1
  while [ $n -le $nmax ]; do
    # initialize loop counter o
    m=1
    while [ $m -le $mmax ]; do
	
  	  # convert indices to hex
	  printf -v memberm "%x" $m
      printf -v membern "%x" $n
	  printf -v membero "%x" $o

      # set index and userconfig name
      RESTART=$RESTARTID
      EXPERIMENT=$ENSEMBLEID"."$membero$membern$memberm
  
      # submit
      echo $EXPERIMENT" / "$RESTART
      subcookie.sh $BASECONFIG $USERCONFIGPATH $EXPERIMENT $YEARS $RESTARTID

      # end loops
      let m=$m+1
    done
    let n=$n+1
  done
  let o=$o+1
done
    