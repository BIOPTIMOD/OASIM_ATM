#!/usr/bin/env bash

# Release build with the Intel compiler inside an OPA chain: cmake and the
# compiler (Intel oneAPI ifx) are taken from the spack environment of the chain,
# in $OPA_HOME/HOST/$OPA_HOSTNAME/spack_env; no module has to be loaded.

DEST_DIR=builds/release_intel

if [ -z "$OPA_HOME" ] || [ -z "$OPA_HOSTNAME" ] ; then
    echo "OPA_HOME and OPA_HOSTNAME must be defined: this script can be used only inside an OPA chain" >&2
    exit 1
fi
SPACK_VIEW=$OPA_HOME/HOST/$OPA_HOSTNAME/spack_env/.spack-env/view
export PATH=$SPACK_VIEW/oneapi-compilers/compiler/latest/bin:$SPACK_VIEW/bin:$PATH

CURR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

rm -rf ${CURR_DIR}/${DEST_DIR}
mkdir -p ${CURR_DIR}/${DEST_DIR}

cd ${CURR_DIR}/${DEST_DIR} && cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_Fortran_COMPILER=$(command -v ifx) ../..
