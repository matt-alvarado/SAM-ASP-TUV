# syntax=docker/dockerfile:1
# SAM-ASP-TUV in a container: builds the model (Intel ifx, Open MPI, NetCDF)
# and runs one fire case (docker/run_case.py). See README, "Running in a
# container".
#
# Build (the common_ASP submodule must be checked out):
#   git submodule update --init
#   docker build -t sam-asp-tuv .
# Options:
#   --build-arg NSUBDOMAINS_X=8   MPI ranks (subdomains in x; must divide 240)
#   --secret id=ca,src=ca.crt     extra CA certificate, for building behind a
#                                 TLS-intercepting proxy (used during the build only)

ARG UBUNTU=ubuntu:24.04

# ---------------------------------------------------------------- toolchain
FROM ${UBUNTU} AS build
ARG DEBIAN_FRONTEND=noninteractive
RUN --mount=type=secret,id=ca,required=false \
    apt-get update -qq && apt-get install -y -qq --no-install-recommends ca-certificates > /dev/null && \
    if [ -s /run/secrets/ca ]; then cp /run/secrets/ca /usr/local/share/ca-certificates/build-ca.crt && update-ca-certificates > /dev/null; fi && \
    sed -i 's/^Types: deb$/Types: deb deb-src/' /etc/apt/sources.list.d/ubuntu.sources && \
    apt-get update -qq && apt-get install -y -qq --no-install-recommends \
        curl bzip2 xz-utils make perl csh dpkg-dev gfortran \
        openmpi-bin libopenmpi-dev libnetcdf-dev libhdf5-dev > /dev/null && \
    rm -rf /var/lib/apt/lists/*
# Intel ifx and icx (conda-forge)
RUN mkdir -p /opt/micromamba && cd /opt/micromamba && \
    curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2 | tar xj bin/micromamba && \
    MAMBA_ROOT_PREFIX=/opt/micromamba ./bin/micromamba create -y -q -p /opt/intel-ifx -c conda-forge \
        ifx_linux-64 dpcpp_linux-64 && \
    rm -rf /opt/micromamba/pkgs
# NetCDF-Fortran rebuilt with ifx (Ubuntu's is built with gfortran)
RUN apt-get update -qq && B=$(mktemp -d) && cd $B && \
    apt-get source -qq --download-only netcdf-fortran && tar xJf netcdf-fortran_*.orig.tar.xz && \
    cd netcdf-fortran-*/ && \
    PATH=/opt/intel-ifx/bin:$PATH FC=ifx F77=ifx CC=icx CPPFLAGS=-I/usr/include \
        LDFLAGS=-L/usr/lib/x86_64-linux-gnu ./configure --prefix=/opt/netcdf-ifx --disable-static -q > /dev/null && \
    cd fortran && PATH=/opt/intel-ifx/bin:$PATH make -j"$(nproc)" -s > /dev/null 2>&1 && make install -s > /dev/null 2>&1 && \
    cd / && rm -rf $B /var/lib/apt/lists/*
RUN printf '%s\n' \
    'export PATH=/opt/intel-ifx/bin:$PATH' \
    'export LIB_SYS=/usr/lib/x86_64-linux-gnu' \
    'export INC_MPI=$LIB_SYS/openmpi/include LIB_MPI_DIR=$LIB_SYS/openmpi/lib' \
    'export INC_NETCDF=/opt/netcdf-ifx/include LIB_NETCDF_DIR=/opt/netcdf-ifx/lib' \
    'export LDFLAGS_TOOLCHAIN="-Wl,-rpath,/opt/intel-ifx/lib -Wl,-rpath-link,/opt/intel-ifx/lib -Wl,-rpath-link,$LIB_SYS -Wl,-rpath,/opt/netcdf-ifx/lib"' \
    > /opt/samasp-env.sh

# ---------------------------------------------------------------- model
WORKDIR /opt/sam-asp-tuv
COPY . .
ARG NSUBDOMAINS_X=16
RUN test -f SRC_TUV/ASP/src/StepASP4SAM.f90 || \
        { echo 'ASP source missing: run "git submodule update --init" before docker build'; exit 1; } && \
    nx=$(sed -n 's/.*nx_gl *= *\([0-9]*\).*/\1/p' SRC_TUV/domain.f90 | head -1) && \
    [ $((nx % NSUBDOMAINS_X)) -eq 0 ] || { echo "NSUBDOMAINS_X=$NSUBDOMAINS_X must divide nx_gl=$nx"; exit 1; } && \
    sed -i "s/^\( *integer, parameter :: nsubdomains_x *= *\)[0-9]*/\1${NSUBDOMAINS_X}/" SRC_TUV/domain.f90 && \
    ny=$(sed -n 's/.*nsubdomains_y *= *\([0-9]*\).*/\1/p' SRC_TUV/domain.f90 | head -1) && \
    echo $((NSUBDOMAINS_X * ny)) > NRANKS && \
    SAMASP_ENV=/opt/samasp-env.sh ./build_ifx.sh > /tmp/build.out 2>&1 || \
        { tail -30 /tmp/build.out; tail -40 _build/sam_build.log 2>/dev/null; exit 1; } && \
    rm -rf _build OBJ && \
    mkdir -p /runtime-libs && \
    ldd SAM_RAD_CAM_MICRO_M2005_w_aero | awk '$3 ~ /^\/opt\// {print $3}' | \
        xargs -I{} cp --parents -L {} /runtime-libs/

# ---------------------------------------------------------------- runtime
FROM ${UBUNTU}
ARG DEBIAN_FRONTEND=noninteractive
RUN --mount=type=secret,id=ca,required=false \
    apt-get update -qq && apt-get install -y -qq --no-install-recommends ca-certificates > /dev/null && \
    if [ -s /run/secrets/ca ]; then cp /run/secrets/ca /usr/local/share/ca-certificates/build-ca.crt && update-ca-certificates > /dev/null; fi && \
    apt-get install -y -qq --no-install-recommends \
        openmpi-bin libnetcdf19t64 python3 python3-numpy python3-netcdf4 python3-pip > /dev/null && \
    pip3 install -q --no-cache-dir --break-system-packages cdsapi && \
    apt-get purge -y -qq python3-pip > /dev/null && apt-get autoremove -y -qq > /dev/null && \
    rm -f /usr/local/share/ca-certificates/build-ca.crt && update-ca-certificates --fresh > /dev/null && \
    rm -rf /var/lib/apt/lists/*
COPY --from=build /runtime-libs/opt /opt
COPY --from=build /opt/sam-asp-tuv /opt/sam-asp-tuv
ENV SAMASP_HOME=/opt/sam-asp-tuv \
    CDSAPI_URL=https://cds.climate.copernicus.eu/api \
    OMPI_MCA_btl_vader_single_copy_mechanism=none
WORKDIR /opt/sam-asp-tuv
VOLUME ["/output"]
ENTRYPOINT ["python3", "/opt/sam-asp-tuv/docker/run_case.py"]
