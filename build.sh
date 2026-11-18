#!/bin/bash
# Build ASP, TUV and SAM with Open MPI and NetCDF, using either the Intel
# oneAPI compilers (ifx/icx) or the GNU compilers (gfortran/gcc, e.g. for
# ARM64 / Apple silicon). ASP and TUV are built out of tree in _build/, then
# SAM with ./Build; the executable is SAM_*_w_aero in this directory.
#
# Usage: ./build.sh [ifx|gnu] [clean]
#   ifx|gnu   compilers (default: $SAMASP_FC, else ifx)
#   clean     remove _build/ and OBJ/ first
# The toolchain paths come from the environment file $SAMASP_ENV (default
# /opt/samasp-env.sh, written by the cloud set-up script docs/cloud_setup.sh
# or by the Dockerfiles); see README, "Build variables".
set -e
cd "$(dirname "$0")"
TOP=$PWD

case $1 in
  ifx) SAMASP_FC=ifx; shift ;;
  gnu|gfortran) SAMASP_FC=gfortran; shift ;;
esac
. "$TOP/toolchain.sh"

if [ "$1" = clean ]; then rm -rf "$TOP/_build" "$TOP/OBJ"; fi

# ASP source is the common_ASP submodule
[ -f SRC_TUV/ASP/src/StepASP4SAM.f90 ] || {
  echo "ASP source missing: run 'git submodule update --init' first"; exit 1; }

mkdir -p _build/asp/ObjectFiles _build/asp/ModuleFiles _build/tuv
echo "== ASP"
(cd _build/asp && make -f "$TOP/SRC_TUV/ASP/Makefile" all \
     ASP_BUILD_DIR="$PWD" ASP_FC=$FC > build.log 2>&1) \
  || { tail -20 _build/asp/build.log; exit 1; }
echo "== TUV"
(cd _build/tuv && make -f "$TOP/SRC_TUV/TUV/Makefile_SAM" all \
     TUV_BUILD_DIR="$PWD" TUV_FC=$FC > build.log 2>&1) \
  || { tail -20 _build/tuv/build.log; exit 1; }

echo "== SAM"
export SAM_FC=$FC SAM_CC=$CC
export ASP_LIBRARY_DIR=$TOP/_build/asp TUV_LIBRARY_DIR=$TOP/_build/tuv
export LIB_MPI="$LIB_MPI_DIR -lmpi_mpifh -lmpi"
export LIB_NETCDF="$LIB_NETCDF_DIR -lnetcdff -L$LIB_SYS -lnetcdf"
export HDF5_HL="$LIB_SYS/hdf5/serial -lhdf5_hl -lhdf5"
export Z="$LIB_SYS -lz -lm"
export SAM_LDFLAGS_EXTRA="$LDFLAGS_TOOLCHAIN"
./Build > _build/sam_build.log 2>&1 || { grep -n "error" _build/sam_build.log | head -20; exit 1; }
ls -l SAM_*_w_aero
