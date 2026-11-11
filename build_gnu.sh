#!/bin/bash
# Build ASP, TUV and SAM with the GNU compilers (gfortran, gcc), Open MPI and
# NetCDF: the build for ARM64 (Apple silicon) and any machine without Intel
# compilers. Same steps and environment file as build_ifx.sh.
#   SAMASP_ENV=/path/to/env.sh ./build_gnu.sh [clean]
SAMASP_FC=gfortran exec "$(dirname "$0")/build_ifx.sh" "$@"
