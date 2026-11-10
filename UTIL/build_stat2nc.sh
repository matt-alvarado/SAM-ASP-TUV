#!/bin/bash
# Build UTIL/stat2nc with the model's toolchain (Intel ifx + NetCDF-Fortran).
# Uses the same environment file as build_ifx.sh:
#   SAMASP_ENV=/opt/samasp-env.sh UTIL/build_stat2nc.sh
set -e
. "${SAMASP_ENV:-/opt/samasp-env.sh}"
cd "$(dirname "$0")/SRC"
ifx -O2 -fixed -extend-source -I. -I"$INC_NETCDF" -o ../stat2nc \
    stat2nc.f hbuf_lib.f cape.f cin.f \
    -L"$LIB_NETCDF_DIR" -L"$LIB_SYS" -lnetcdff -lnetcdf $LDFLAGS_TOOLCHAIN
echo "built $(cd .. && pwd)/stat2nc"
