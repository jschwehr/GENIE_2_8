#!/bin/bash
echo "Setting GENIE environment variables..."
export GENIEBASE=/home/jackie/work/genie
export GENIE=$GENIEBASE/GENIE_2_8
export PYTHIA6=/software/GENIE/support/pythia6/v6_424/lib
export ROOTSYS=/software/GENIE/support/root.5.34.14
export LOG4CPP_INC=/usr/local/include/log4cpp
export LOG4CPP_LIB=/usr/local/lib
export LHAPATH=/software/GENIE/support/lhapdf
export LHAPDF_INC=/software/GENIE/support/lhapdf/include
export LHAPDF_LIB=/software/GENIE/support/lhapdf/lib
export XSECSPLINEDIR=$GENIEBASE/data
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LHAPDF_LIB:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LOG4CPP_LIB:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYTHIA6:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH
export PATH=$GENIE/bin:$ROOTSYS/bin:$PATH
unset GENIEBASE
