#!/bin/bash

. /etc/os-release

#please modify your configuration here
R_PACKAGES_PATH=/data_featherstone/christiane/Rpackages     #the path to store the packages
export R_VERSION=4.2.1

####setup####
if [ ! -d /sw/R/${VERSION_CODENAME}/${R_VERSION} ]; then
    echo "Error: Please check your R version"
    echo "Installed versions:"
    ls /sw/R/${VERSION_CODENAME}/
    exit
fi
export PATH=/sw/R/${VERSION_CODENAME}/${R_VERSION}/bin/:$PATH
export R_LIBS_USER=${R_PACKAGES_PATH}/${R_VERSION}
export R_HOME=/sw/R/${VERSION_CODENAME}/${R_VERSION}

if [ ! -d ${R_PACKAGES_PATH}/${R_VERSION} ]; then
    mkdir -p ${R_PACKAGES_PATH}/${R_VERSION} || exit
fi

if [ -z $LD_LIBRARY ]; then
    export LD_LIBRARY=/sw/R/${VERSION_CODENAME}/${R_VERSION}/lib/R/lib/
else
    export LD_LIBRARY=/sw/R/${VERSION_CODENAME}/${R_VERSION}/lib/R/lib/:$LD_LIBRARY
fi    

export OMP_NUM_THREADS=1 # number of threads to be used while running the job

### input arguments to be used as input
rown=$1

Rscript /xx/xxx/SC/CAM_TUM/eqtl/Targeted/Scripts/Permutation.R $rown
