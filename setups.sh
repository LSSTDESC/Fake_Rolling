#!/bin/bash

##source /sps/lsst/Library/stack_v11_0/loadLSST.bash
#export LSSTSW=/sps/lsst/Library/lsstsw
#export EUPS_PATH=$LSSTSW/stack
#source $LSSTSW/bin/setup.sh
##setup pipe_tasks
#setup sconsUtils 10.0
#C++11

export PATH=/opt/rh/devtoolset-3/root/usr/bin:${PATH}

#export LSSTSW=/sps/lsst/Library/lsstsw
export LSSTSW=/sps/lsst/Library/new/lsstsw
export EUPS_PATH=$LSSTSW/stack
source $LSSTSW/bin/setup.sh

dir_tag=/sps/lsst/data/dev/pgris/sims_operations

export PYTHONPATH=${dir_tag}/lib/python2.7/site-packages:$PYTHONPATH

setup -k -r ${dir_tag}/pykg_config
#setup -k -r ${dir_tag}/healpy
setup -k -r ${dir_tag}/python_future
setup -k -r ${dir_tag}/sims_data
setup -k -r ${dir_tag}/sims_utils
#setup -k -r ${dir_tag}/sims_catalogs_generation
setup -k -r ${dir_tag}/sims_sed_library
setup -k -r ${dir_tag}/sims_dustmaps
setup -k -r ${dir_tag}/throughputs
setup -k -r ${dir_tag}/sims_maps
setup -k -r ${dir_tag}/sims_photUtils
setup -k -r ${dir_tag}/sims_coordUtils
setup -k -r ${dir_tag}/sqlalchemy
setup -k -r ${dir_tag}/sims_catalogs
setup -k -r ${dir_tag}/sims_maf
export PYTHONPATH=/sps/lsst/data/dev/pgris/sims_operations/sims_operations/python/lsst/sims/operations:$PYTHONPATH
#setup -k -r ../sims_operations