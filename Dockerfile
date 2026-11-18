# syntax=docker/dockerfile:1
# SAM-ASP-TUV in a container: builds the model (Intel ifx, Open MPI, NetCDF)
# and runs one fire case (docker/run_case.py). See README, "Running in a
# container".
#
# Build (the common_ASP submodule must be checked out):
#   git submodule update --init
#   docker build -t sam-asp-tuv .
# Options:
#   --build-arg NSUBDOMAINS_X=8   MPI ranks (subdomains in x; must divide nx_gl)
#   --build-arg NTRACERS=1486     SAM tracer count for decks other than the default
#                                 (per case in cases/README.md; the model prints
#                                 the value a deck set needs)
#   --secret id=ca,src=ca.crt     extra CA certificate, for building behind a
#                                 TLS-intercepting proxy (used during the build only)

# fully qualified, so Podman/Buildah (RHEL) need no short-name lookup
ARG UBUNTU=docker.io/library/ubuntu:24.04

# ---------------------------------------------------------------- toolchain
FROM ${UBUNTU} AS build
ARG DEBIAN_FRONTEND=noninteractive
# The toolchain is that of the cloud development environment (docs/cloud_setup.sh):
# Open MPI and NetCDF-C from Ubuntu, Intel ifx/icx from conda-forge,
# NetCDF-Fortran rebuilt with ifx (Ubuntu's is built with gfortran), and the
# environment file /opt/samasp-env.sh that build.sh reads.
COPY docker/add_build_ca.sh docs/cloud_setup.sh /tmp/
RUN --mount=type=secret,id=ca,required=false \
    apt-get update -qq && apt-get install -y -qq --no-install-recommends \
        ca-certificates curl make perl dpkg-dev > /dev/null && \
    sh /tmp/add_build_ca.sh && \
    bash /tmp/cloud_setup.sh && \
    rm -rf /opt/micromamba/pkgs /var/lib/apt/lists/*

# ---------------------------------------------------------------- model
WORKDIR /opt/sam-asp-tuv
COPY . .
ARG NSUBDOMAINS_X=16
ARG NTRACERS=
RUN docker/build_model.sh ifx && \
    mkdir -p /runtime-libs && \
    ldd SAM_RAD_CAM_MICRO_M2005_w_aero $(sed -n 's/^\//UTIL\//p' UTIL/.gitignore) | awk '$3 ~ /^\/opt\// {print $3}' | \
        sort -u | xargs -I{} cp --parents -L {} /runtime-libs/

# ---------------------------------------------------------------- runtime
FROM ${UBUNTU}
ARG DEBIAN_FRONTEND=noninteractive
COPY docker/add_build_ca.sh /tmp/
RUN --mount=type=secret,id=ca,required=false \
    apt-get update -qq && apt-get install -y -qq --no-install-recommends ca-certificates > /dev/null && \
    sh /tmp/add_build_ca.sh && \
    apt-get install -y -qq --no-install-recommends \
        openmpi-bin libnetcdf19t64 python3 python3-numpy python3-netcdf4 python3-pip > /dev/null && \
    pip3 install -q --no-cache-dir --break-system-packages cdsapi && \
    apt-get purge -y -qq python3-pip > /dev/null && apt-get autoremove -y -qq > /dev/null && \
    sh /tmp/add_build_ca.sh --remove && rm /tmp/add_build_ca.sh && \
    rm -rf /var/lib/apt/lists/*
COPY --from=build /runtime-libs/opt /opt
COPY --from=build /opt/sam-asp-tuv /opt/sam-asp-tuv
ENV SAMASP_HOME=/opt/sam-asp-tuv \
    CDSAPI_URL=https://cds.climate.copernicus.eu/api \
    OMPI_MCA_btl_vader_single_copy_mechanism=none
WORKDIR /opt/sam-asp-tuv
VOLUME ["/output"]
ENTRYPOINT ["python3", "/opt/sam-asp-tuv/docker/run_case.py"]
