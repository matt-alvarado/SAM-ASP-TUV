#!/bin/bash
# Build the SAM post-processing tools in UTIL/ (stat2nc, bin2D2nc, ...) from
# UTIL/SRC with the model's toolchain (Intel ifx + NetCDF-Fortran). Uses the
# same environment file as build_ifx.sh:
#   SAMASP_ENV=/opt/samasp-env.sh UTIL/build_util.sh            # all tools
#   SAMASP_ENV=/opt/samasp-env.sh UTIL/build_util.sh stat2nc    # some tools
# The binaries are written to UTIL/ (git-ignored; no prebuilt ones are kept).
# SAMASP_FC=gfortran builds them with GNU Fortran (as build_gnu.sh does the model).
set -e
. "${SAMASP_ENV:-/opt/samasp-env.sh}"
case ${SAMASP_FC:-ifx} in
  ifx) FC="ifx -O2 -fixed -extend-source"; EXTRA= ;;
  gfortran) FC="gfortran -O2 -ffixed-form -ffixed-line-length-132 -fallow-argument-mismatch -std=legacy"
            EXTRA=iargc_gnu.f ;;
  *) echo "SAMASP_FC must be ifx or gfortran"; exit 1 ;;
esac
TOOLS="2Dbin2nc 2Dbin2nc_mean 2Dbin2nc_sep 2Dbin_sep2one 2Dcom2nc 2Dcom2nc_mean 2Dcom2nc_sep 2Dcom_sep2one bin2D2nc bin3D2nc bin3D2nc_mean bin3D_sep2one com2D2nc com3D2bin com3D2nc com3D2nc_mean com3D2nc_sep com3D_sep2one glue_movie_raw isccp2nc misr2nc modis2nc stat2nc"
cd "$(dirname "$0")/SRC"
for t in ${@:-$TOOLS}; do
    $FC -I. -I"$INC_NETCDF" -o ../$t \
        $t.f hbuf_lib.f cape.f cin.f $EXTRA \
        -L"$LIB_NETCDF_DIR" -L"$LIB_SYS" -lnetcdff -lnetcdf $LDFLAGS_TOOLCHAIN
    echo "built UTIL/$t"
done
