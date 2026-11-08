#!/bin/bash
# Run the short SAM-ASP-TUV coupling test (RUNS/asptest) with the ifx build.
#
# The run happens in a copy of the repository, because a run overwrites files
# that are committed (SRC_TUV/ASP/OutputFiles/*, CaseName, RUNS/*/*.nml).
#
# Usage: ./build_ifx.sh && ./run_asptest_ifx.sh [run_dir]
#   run_dir defaults to ./_run_asptest (git-ignored).
# Takes about 12 minutes with 16 ranks on 4 cores. A clean finish exits with
# code 9 (SAM's "reached nstop" code, see SRC_TUV/task_util_MPI.f90).
set -e
cd "$(dirname "$0")"
TOP=$PWD
RUN=${1:-$TOP/_run_asptest}
EXE=SAM_RAD_CAM_MICRO_M2005_w_aero
[ -x "$EXE" ] || { echo "Build first: ./build_ifx.sh"; exit 1; }

mkdir -p "$RUN"
tar --exclude=.git --exclude=_build --exclude=OBJ --exclude=_run_asptest \
    -cf - . | tar -xf - -C "$RUN"
cd "$RUN"
echo asptest > CaseName
mkdir -p OUT_2D OUT_3D OUT_MOMENTS OUT_STAT OUT_MOVIES RESTART
rm -f OUT_STAT/asptest_asptest.stat   # SAM refuses to overwrite it

ulimit -s unlimited   # TUV and ASP use large local arrays
set +e
mpirun.openmpi --allow-run-as-root --oversubscribe -np 16 ./$EXE > run.out 2>&1
rc=$?
set -e
echo "exit code $rc (9 = finished normally); log: $RUN/run.out"
# ASP calls are counted from a per-call message printed only with the
# namelist option samasp_verbose = .true.
nasp=$(grep -a -c 'After SAM_wrapper' run.out)
[ "$nasp" -gt 0 ] || nasp="(set samasp_verbose to count)"
echo "ASP calls: $nasp," \
     "last step: $(grep -a 'NSTEP =' run.out | tail -1 | awk '{print $3}')," \
     "Fortran runtime errors: $(grep -a -c forrtl run.out)" \
     "(a 'process killed (SIGTERM)' at the very end is Open MPI stopping the" \
     "other ranks after SAM's normal exit code 9)"
