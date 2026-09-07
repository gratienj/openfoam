#!/bin/bash

# Purge all loaded modules
#module purge

# Update MODULEPATH
#if [ -z "${BASH_SOURCE[0]}" ]; then
#  export MODULEPATH="$(dirname $(dirname `readlink -f $0`))/easybuild/modules/all/Core":$MODULEPATH
#else
#  export MODULEPATH="$(dirname $(dirname `readlink -f ${BASH_SOURCE[0]}`))/easybuild/modules/all/Core":$MODULEPATH
#fi

export MODULEPATH=/soft/irsrvsoft1/expl/eb/r11/toolchains/el/9/x86_64/easybuild/modules/all/Core:$MODULEPATH
# Toolchain
module load foss/2025b
#export LD_LIBRARY_PATH=$PWD/libefa:$LD_LIBRARY_PATH

