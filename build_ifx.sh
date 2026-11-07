#!/bin/bash
# Build ASP, TUV and SAM with the Intel oneAPI ifx/icx compilers, Open MPI
# and NetCDF, as installed by the "SAM-ASP" cloud environment's setup script
# (which writes /opt/samasp-env.sh).
#
# ASP and TUV are built out of tree in _build/, so the prebuilt ifort
# libraries, objects and modules committed under SRC_TUV/ are left untouched.
# The rotor build (./Build with ifort) is unchanged.
#
# Usage: ./build_ifx.sh            build everything
#        ./build_ifx.sh clean      remove _build/ and OBJ/ first
set -e
cd "$(dirname "$0")"
TOP=$PWD

ENV_SCRIPT=${SAMASP_ENV:-/opt/samasp-env.sh}
[ -f "$ENV_SCRIPT" ] || { echo "Missing $ENV_SCRIPT (toolchain setup)"; exit 1; }
source "$ENV_SCRIPT"

if [ "$1" = clean ]; then rm -rf "$TOP/_build" "$TOP/OBJ"; fi

# ASP source is the common_ASP submodule
[ -f SRC_TUV/ASP/src/StepASP4SAM.f90 ] || {
  echo "ASP source missing: run 'git submodule update --init' first"; exit 1; }

mkdir -p _build/asp/ObjectFiles _build/asp/ModuleFiles _build/tuv
echo "== ASP"
(cd _build/asp && make -f "$TOP/SRC_TUV/ASP/Makefile" all \
     ASP_BUILD_DIR="$PWD" ASP_FC=ifx > build.log 2>&1) \
  || { tail -20 _build/asp/build.log; exit 1; }
echo "== TUV"
(cd _build/tuv && make -f "$TOP/SRC_TUV/TUV/Makefile_SAM" all \
     TUV_BUILD_DIR="$PWD" TUV_FC=ifx > build.log 2>&1) \
  || { tail -20 _build/tuv/build.log; exit 1; }

echo "== SAM"
export SAM_FC=ifx SAM_CC=icx
export ASP_LIBRARY_DIR=$TOP/_build/asp TUV_LIBRARY_DIR=$TOP/_build/tuv
export LIB_MPI="$LIB_MPI_DIR -lmpi_mpifh -lmpi"
export LIB_NETCDF="$LIB_NETCDF_DIR -lnetcdff -L$LIB_SYS -lnetcdf"
export HDF5_HL="$LIB_SYS/hdf5/serial -lhdf5_hl -lhdf5"
export Z="$LIB_SYS -lz -lm"
export SAM_LDFLAGS_EXTRA="$LDFLAGS_TOOLCHAIN"
./Build > _build/sam_build.log 2>&1 || { grep -n "error" _build/sam_build.log | head -20; exit 1; }
ls -l SAM_*_w_aero
