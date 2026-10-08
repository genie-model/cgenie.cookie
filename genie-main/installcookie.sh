#!/bin/bash -e
#
#####################################################################
### SCIPT TO COMPLETE INSTALLATION OF CGENIE.COOKIE #################
#####################################################################
#
echo ""
#
# ---------------------------------------------------------------------
# (1) LOOKUP TABLES
# ---------------------------------------------------------------------
#
## unpack lookup tables
#echo " > Unpacking lookup tables ..."
#tar xzf ../genie-sedgem/data/input/lookup_calcite_4.tar.gz -C ../genie-sedgem/data/input/
#tar xzf ../genie-sedgem/data/input/lookup_opal_5.tar.gz -C ../genie-sedgem/data/input/
#
# ---------------------------------------------------------------------
# (2) PYTHON VERSION
# ---------------------------------------------------------------------
#
# check host name ...
# ... and copy python2 into genie-main 
# if that will make the world a better place
SERVER=$(hostname -s)
if [[ "$SERVER" == "sterling" || "$SERVER" == "eevee" ]]; then
   echo '>> copying python2 files for: '$SERVER
   cp src/python2/* .
else
   echo ' > Current host: '$SERVER' is assumed to be python3-friendly'
fi
#
## determine python2 vs.3 and copy apropriate set of files
#if [[ -z "$(which python3 2>/dev/null)" ]]
#  cp python3/* .
#else
#  cp python3/* .
#fi
#
# ---------------------------------------------------------------------
# END
# ---------------------------------------------------------------------
# 
echo ""
echo ">> cookie installation complete <<"
echo ""
#
#####################################################################
#####################################################################
