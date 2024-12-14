#!/bin/bash -e
#
########################################################################
### SCIPT TO COMPLETE INSTALLATION OF CGENIE.COOKIE ####################
########################################################################
#
echo ""
#
########################################################################
# determine python version ... copy appropriate file versions



#
########################################################################
# unpack lookup tables
echo ">> Unpacking lookup tables ..."
tar xzf ../genie-sedgem/data/input/lookup_calcite_4.tar.gz -C ../genie-sedgem/data/input/
tar xzf ../genie-sedgem/data/input/lookup_opal_5.tar.gz -C ../genie-sedgem/data/input/
# 
########################################################################
# END
echo "   *** muffin installation complete ***"
echo ""
#
