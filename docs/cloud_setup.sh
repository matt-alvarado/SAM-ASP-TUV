#!/bin/bash
# SAM-ASP-TUV toolchain: Intel ifx (conda-forge), Open MPI, NetCDF (C from apt, Fortran rebuilt with ifx)
#
# Copy of the Claude Code cloud environment's setup script (paste it into the
# environment's settings: cloud environment menu > Edit > Setup script). It
# writes /opt/samasp-env.sh, which build_ifx.sh / build_gnu.sh source.
# csh is needed by SAM's Build script.
#
# Network access for that environment: "Limited" with package managers
# allowed, plus these hosts:
#   object-store.os-api.cci2.ecmwf.int   ERA5 downloads (CDS result files)
#   registry.access.redhat.com, cdn-ubi.redhat.com, cdn01.quay.io
#                                        UBI 9 images and packages
#   zenodo.org                           Zenodo deck archives
# The Copernicus CDS key is set as a network secret for
# cds.climate.copernicus.eu (never put it in this file or the repository).
set -e
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq openmpi-bin libopenmpi-dev libnetcdf-dev netcdf-bin libhdf5-dev gfortran xz-utils bzip2 csh >/dev/null
if [ ! -x /opt/intel-ifx/bin/ifx ]; then
  mkdir -p /opt/micromamba && cd /opt/micromamba
  curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2 | tar xj bin/micromamba
  MAMBA_ROOT_PREFIX=/opt/micromamba ./bin/micromamba create -y -q -p /opt/intel-ifx -c conda-forge ifx_linux-64 dpcpp_linux-64
  cd - >/dev/null
fi
if [ ! -f /opt/netcdf-ifx/lib/libnetcdff.so ]; then
  B=$(mktemp -d) && cd $B
  sed -i 's/^Types: deb$/Types: deb deb-src/' /etc/apt/sources.list.d/ubuntu.sources && apt-get update -qq
  apt-get source -qq --download-only netcdf-fortran && tar xJf netcdf-fortran_*.orig.tar.xz && cd netcdf-fortran-*/
  PATH=/opt/intel-ifx/bin:$PATH FC=ifx F77=ifx CC=icx CPPFLAGS=-I/usr/include LDFLAGS=-L/usr/lib/x86_64-linux-gnu ./configure --prefix=/opt/netcdf-ifx --disable-static -q >/dev/null
  (cd fortran && PATH=/opt/intel-ifx/bin:$PATH make -j4 -s >/dev/null 2>&1 && make install -s >/dev/null 2>&1)
  cd / && rm -rf $B
fi
cat > /opt/samasp-env.sh <<'ENV'
export PATH=/opt/intel-ifx/bin:$PATH
export LIB_SYS=/usr/lib/x86_64-linux-gnu
export INC_MPI=$LIB_SYS/openmpi/include LIB_MPI_DIR=$LIB_SYS/openmpi/lib
export INC_NETCDF=/opt/netcdf-ifx/include LIB_NETCDF_DIR=/opt/netcdf-ifx/lib
export LDFLAGS_TOOLCHAIN="-Wl,-rpath,/opt/intel-ifx/lib -Wl,-rpath-link,/opt/intel-ifx/lib -Wl,-rpath-link,$LIB_SYS -Wl,-rpath,/opt/netcdf-ifx/lib"
export LIBS_TOOLCHAIN="-L$LIB_NETCDF_DIR -lnetcdff -L$LIB_SYS -lnetcdf -L$LIB_MPI_DIR -lmpi_mpifh -lmpi"
alias mpirun='mpirun.openmpi --allow-run-as-root --oversubscribe'
ENV
