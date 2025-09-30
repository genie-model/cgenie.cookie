#!/bin/bash -e
#
#####################################################################
### SCIPT TO COMPLETE INSTALLATION OF CGENIE.COOKIE #################
#####################################################################
#
echo ""
#
# unpack lookup tables
echo ">> Unpacking lookup tables ..."
tar xzf ../genie-sedgem/data/input/lookup_calcite_4.tar.gz -C ../genie-sedgem/data/input/
tar xzf ../genie-sedgem/data/input/lookup_opal_5.tar.gz -C ../genie-sedgem/data/input/
# 
echo "   *** cookie installation complete ***"
echo ""
#

# determine python2 vs.3 and copy apropriate set of files
if [[ -z "$(which python3 2>/dev/null)" ]]
  cp python3/* .
else
  cp python3/* .
fi
#
