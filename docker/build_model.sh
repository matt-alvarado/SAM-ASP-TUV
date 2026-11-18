#!/bin/bash
# Model stage of the Dockerfiles, run in the repository directory:
#   NSUBDOMAINS_X=16 NTRACERS= docker/build_model.sh ifx|gnu
# 1. sets nsubdomains_x (and ntracers, if NTRACERS is set) in
#    SRC_TUV/domain.f90; NSUBDOMAINS_X must divide nx_gl there;
# 2. writes NRANKS (MPI ranks = nsubdomains_x * nsubdomains_y), read by
#    docker/run_case.py;
# 3. builds the model (build.sh) and the UTIL tools with that compiler and the
#    environment file $SAMASP_ENV (default /opt/samasp-env.sh), then removes
#    the build tree. On failure the end of the build logs is printed.
set -e
cd "$(dirname "$0")/.."
[ $# -eq 1 ] || { echo "usage: $0 ifx|gnu" >&2; exit 1; }
test -f SRC_TUV/ASP/src/StepASP4SAM.f90 || {
    echo 'ASP source missing: run "git submodule update --init" before the build' >&2; exit 1; }

DOMAIN=SRC_TUV/domain.f90
param() { sed -n "s/^ *integer, parameter :: $1 *= *\([0-9]*\).*/\1/p" $DOMAIN | head -1; }
set_param() { sed -i "s/^\( *integer, parameter :: $1 *= *\)[0-9]*/\1$2/" $DOMAIN; }

NSUBDOMAINS_X=${NSUBDOMAINS_X:?set NSUBDOMAINS_X}
nx=$(param nx_gl)
[ $((nx % NSUBDOMAINS_X)) -eq 0 ] || {
    echo "NSUBDOMAINS_X=$NSUBDOMAINS_X must divide nx_gl=$nx" >&2; exit 1; }
set_param nsubdomains_x "$NSUBDOMAINS_X"
[ -z "$NTRACERS" ] || set_param ntracers "$NTRACERS"
echo $((NSUBDOMAINS_X * $(param nsubdomains_y))) > NRANKS

case $1 in
  ifx) export SAMASP_FC=ifx ;;
  gnu) export SAMASP_FC=gfortran ;;
  *) echo "usage: $0 ifx|gnu" >&2; exit 1 ;;
esac
./build.sh > /tmp/build.out 2>&1 || {
    tail -30 /tmp/build.out; tail -40 _build/sam_build.log 2>/dev/null; exit 1; }
UTIL/build_util.sh > /dev/null
rm -rf _build OBJ
