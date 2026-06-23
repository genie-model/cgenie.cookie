#!/bin/bash -e
#
#####################################################################
### SCIPT TO META-CONFIGURE AND SUBMIT CGENIE.cookie ################
#####################################################################
#
echo ""
# ---------------------------------------------------------------------
# (0) USER OPTIONS
# ---------------------------------------------------------------------
#####################################################################
# CHANGE THIS FOR SPECIFIC QUEUES
# specific any particular queue to be used (empty string otherwise)
QUEUE='-q dog.q'
# CHANGE THIS FOR INSTALLATIONS OTHER THAN IN $HOME
# SET THE SAME AS IN user.mak AND user.sh
# set home directory
HOMEDIR=$HOME
# set output directory
OUTPUTDIR=$HOMEDIR/cgenie_output
#####################################################################
# ---------------------------------------------------------------------
# (1) GET PASSED PARAMETERS
# ---------------------------------------------------------------------
# [1] base configuration ID
if [ -z "$1" ]; then
    echo "Usage: '$1' 1st parameter must be the base configuration ID"
    exit 65
  else
    MODELID="$1"
    SUB1=$1
fi
# [2] set user (experiment) configuration file directory
if [ -z "$2" ]; then
    echo "Usage: '$2' 2nd parameter must be the user (experiment) configuration file directory"
    exit 65
  else
    GOINDIR=$HOMEDIR/cgenie.cookie/genie-userconfigs/"$2"
    SUB2=$2
fi
# [3] set run ID (input run ID (= user configuration file name))
if [ -z "$3" ]; then
    echo "Usage: '$3' 3rd parameter must be the run ID"
    exit 65
  else
    RUNID="$3"
    SUB3=$3
fi
if [ $(expr length "$3") -gt 95 ] ; then
    echo "Usage: '$3' 3rd parameter must be less than 96 characters in length"
    exit 65
fi
GOIN=$GOINDIR/$RUNID
# [4] set run duration
if [ -z "$4" ]; then
    echo "Usage: '$4' 4th parameter must be the run length (years)"
    exit 65
  else
    RUNLENGTH="$4"
    SUB4=$4
fi
# [5] restart path (optional)
if [ -n "$5" ]; then
  RESTARTPATH=$OUTPUTDIR/"$5"
    if [ $(expr length "$5") -gt 95 ] ; then
        echo "Usage: '$5' 5th parameter must be less than 96 characters in length"
        exit 65
    fi
    SUB5=$5
else
    SUB5=""
fi
# ---------------------------------------------------------------------
# (2) SET LOCAL FILE AND DIRECTORY NAMES
# ---------------------------------------------------------------------
#
OUTPUTPATH=$OUTPUTDIR/$RUNID
CONFIGPATH=$HOMEDIR/cgenie.cookie/genie-baseconfigs
CONFIGNAME=$RUNID".config"
BINARYPATH=$HOMEDIR/cgenie.cookie/genie-main
RESTARTNAME="rst.1"
# ---------------------------------------------------------------------
# (3) CHECK PARAMETERS
# ---------------------------------------------------------------------
#
echo ">> Checking parameters ..."
echo ""
#
# NOTE: deal with ".config" being accidently included in the run command
#
CONFIGEXT=${MODELID: -7}
if [[ "$CONFIGEXT" == ".config" ]]; then
    echo "   #0 Removing .config from base configuration name (before adding it back again later ...): "
    echo $MODELID
    MODELID=${MODELID:: -7}
fi
if test -e $CONFIGPATH/$MODELID".config"
then
    echo "   #1 Experiment base configuration: "
    echo $CONFIGPATH/$MODELID".config"
    echo " found :)"
    SUB1=$MODELID".config"
else
    echo "   #1 Experiment base configuration: "
    echo $CONFIGPATH/$MODELID".config"
    echo " CANNOT be found :("
    exit 1
fi
if test -d $GOINDIR
then
    echo "   #2 Experiment user configuration file directory: "
    echo $GOINDIR
    echo " found :)"
else
    echo "   #2 Experiment user configuration file directory: "
    echo $GOINDIR
    echo " CANNOT be found :("
    exit 1
fi
if test -e $GOIN
then
    echo "   #3 User-config file: "
    echo $GOIN
    echo " found :)"
else
    echo "   #3 User-config file: "
    echo $GOIN
    echo " CANNOT be found :("
    exit 1
fi
if [ $(expr length "$3") -gt 127 ] ; then
    echo "Usage: '$3' 3rd parameter must be less than 128 characters in length"
    exit 65
fi
if [ -n "$5" ]; then
    if test -e $RESTARTPATH
    then
        echo "   #5 Restart experiment: "
        echo $RESTARTPATH
        echo " found :)"
    else
        echo "   #5 Restart experiment: "
        echo $RESTARTPATH
        echo " CANNOT be found :("
        exit 1
    fi
    if [ $(expr length "$5") -gt 127 ] ; then
        echo "Usage: '$5' 5th parameter must be less than 128 characters in length"
        exit 65
    fi
else
    echo "   #5 NO restart specified"
fi
cp -f $CONFIGPATH/$MODELID".config" $CONFIGPATH/$CONFIGNAME
# ---------------------------------------------------------------------
# (4) SET MODEL TIME-STEPPING
# ---------------------------------------------------------------------
# extract ocean (lon,lat) dimension
LONS=$(grep -o '$(DEFINE)GOLDSTEINNLONS=..\>' $CONFIGPATH/$MODELID".config" | sed -e s/.*=//)
LATS=$(grep -o '$(DEFINE)GOLDSTEINNLATS=..\>' $CONFIGPATH/$MODELID".config" | sed -e s/.*=//)
# extract ocean depth dimentsion
# NOTE: test for empty as a single digit ... then try 2
# NOTE: if the first digit is a zero ... remove it ...
LEVS=$(grep -o '$(DEFINE)GOLDSTEINNLEVS=..\>' $CONFIGPATH/$MODELID".config" | sed -e s/.*=//)
if [ -z "$LEVS" ]; then
    LEVS=$(grep -o '$(DEFINE)GOLDSTEINNLEVS=.\>' $CONFIGPATH/$MODELID".config" | sed -e s/.*=//)
else
    if [[ ${LEVS:0:1} == "0" ]]; then LEVS=${LEVS:1:2}; fi
fi
# report
# echo "   Requested model resolution: "$LONS"x"$LATS"x"$LEVS
# ---------------------------------
# (9) GO!
# ---------------------------------
# extract number of tracers -- generalized by Gemini
#TRACERS=$(grep -o '$(DEFINE)GOLDSTEINNTRACS=..\>' $CONFIGPATH/$MODELID".config" | sed -e s/.*=//)
TRACERS=$(grep -v '^[[:space:]]*#' $CONFIGPATH/$MODELID".config" | grep -o '$(DEFINE)GOLDSTEINNTRACS=[0-9]\+' | sed -e 's/.*=//' | tr -d '\r' | head -n 1)
# create array size #
ndim="$(echo "$LONS*$LATS*$LEVS*$TRACERS" | bc -l)"
# test for change of base-config
if test -e 'cookie.current_grid_size.txt' 
then
    cookie_n=$(<'cookie.current_grid_size.txt')
    if [ "$cookie_n" != "$ndim" ]; then
        echo ""
        echo ">> Your configuration requires a different array size to that used previously, so ... "
        echo "   ... we are going to remove the model executable and re-compile from the FORTRAN source code."
        echo ""
        sleep 5
        make cleanall
        ./genie_example.job -O -f $CONFIGPATH/$CONFIGNAME -x
        echo "$ndim" > 'cookie.current_grid_size.txt'
    fi
else
    echo ""
    echo ">> No record of last ocean array size used (file: cookie.current_grid_size.txt) ... COMPILING ..."
    echo ""
    sleep 5
    make cleanall
    ./genie_example.job -O -f $CONFIGPATH/$CONFIGNAME -x
    echo "$ndim" > 'cookie.current_grid_size.txt'
fi
# Submit job ...
echo ""
qsub $QUEUE -j y -o cgenie_log -V -S /bin/bash runcookie.sh $SUB1 $SUB2 $SUB3 $SUB4 $SUB5
echo ""
echo ">> ... waiting for it to appear on the queue (I am just guessing how long that will be and assuming 15s) ... "
echo ""
sleep 15
qstat -f
echo ""
# ---------------------------------
# ---------------------------------
#
