# Sourced by build.sh and UTIL/build_util.sh: loads the toolchain environment
# file ($SAMASP_ENV, default /opt/samasp-env.sh) and sets FC and CC from
# $SAMASP_FC (ifx, the default, or gfortran; "gnu" means gfortran).
ENV_SCRIPT=${SAMASP_ENV:-/opt/samasp-env.sh}
[ -f "$ENV_SCRIPT" ] || { echo "Missing $ENV_SCRIPT (toolchain set-up; see README, Build variables)" >&2; exit 1; }
. "$ENV_SCRIPT"
case ${SAMASP_FC:-ifx} in
  ifx) FC=ifx CC=icx ;;
  gfortran|gnu) FC=gfortran CC=gcc ;;
  *) echo "SAMASP_FC must be ifx or gfortran" >&2; exit 1 ;;
esac
