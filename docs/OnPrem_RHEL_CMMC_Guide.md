# SAM-ASP-TUV on an on-prem RHEL VM in a CMMC Level 2 enclave

This guide covers installing and running the SAM-ASP-TUV container on an
on-premises Linux virtual machine (VM) inside a CMMC Level 2 enclave, that
is, a system that handles Controlled Unclassified Information (CUI) and must
meet the 110 requirements of NIST SP 800-171 rev 2.

It uses the **RHEL 9 container image** (`Dockerfile.ubi9`, built on Red
Hat's Universal Base Image 9), so the container runs the same operating
system as the host. The model, its settings and its outputs are the same as
with the Ubuntu-based image described in the
[container user's guide](Container_User_Guide.md), which explains the model
settings, ASP decks, ERA5, fire timing and outputs. Read that guide for
everything about the model itself; where it says `docker`, use `podman`.
For an Apple silicon Mac, see the [Mac guide](Mac_AppleSilicon_CMMC_Guide.md).

> **This is not a compliance determination.** CMMC sets requirements, not
> products. What is allowed on your VM is decided by your organization's
> System Security Plan (SSP) and its security officer (ISSO/ISSM). The
> control numbers below (e.g. 3.4.8) point to the NIST SP 800-171
> requirements each step touches, so you and the ISSO can find them. Where
> your SSP says otherwise, follow your SSP.

Contents:
1. [Assumptions](#1-assumptions)
2. [What to agree with your ISSO first](#2-what-to-agree-with-your-isso-first)
3. [Machine size: CPU, memory, GPU](#3-machine-size-cpu-memory-gpu)
4. [Administrator setup (once)](#4-administrator-setup-once)
5. [Getting the image](#5-getting-the-image)
6. [User setup (once per user)](#6-user-setup-once-per-user)
7. [Running a case](#7-running-a-case)
8. [Handling inputs and outputs](#8-handling-inputs-and-outputs)
9. [Records, updates and removal](#9-records-updates-and-removal)
10. [Troubleshooting](#10-troubleshooting)
11. [What was tested](#11-what-was-tested)

---

## 1. Assumptions

This guide is written for the most common enclave setup. If yours differs,
the sections note what changes.

| Item | Assumed |
|---|---|
| VM | Inside the CUI enclave; 16 vCPUs (AMD), 64 GB RAM; NVIDIA L4 GPU with CUDA 13.4 (not used by the model, see 3.3) |
| Host OS | Red Hat Enterprise Linux 9, hardened (DISA STIG or CIS), **FIPS mode on**, SELinux enforcing |
| Container OS | **RHEL 9 as well**: the image is built on Red Hat's Universal Base Image 9 (UBI 9), or on your organization's hardened UBI 9 image (`Dockerfile.ubi9`) |
| Container runtime | **Podman, rootless** (no Docker daemon, no root for the person running the model) |
| Privileges | An administrator installs software; model users have no general `sudo` (3.1.5, 3.1.6) |
| Network | Outbound internet through the enclave's proxy, to allowed hosts only (3.13.1, 3.13.6); the proxy may inspect TLS with an organization CA certificate |
| Data | ERA5 meteorology and the model code are public. ASP decks and model outputs **may** be CUI, depending on your project; this guide treats them as if they are |

---

## 2. What to agree with your ISSO first

Bring this list to your ISSO before installing. Each item needs a decision,
and usually an entry in the SSP or the software/configuration records.

| # | Item | Why | Requirement |
|---|---|---|---|
| 1 | **Approve the software**: Podman (`container-tools`) on the VM, the base image (Red Hat UBI 9, or your hardened UBI 9), and the SAM-ASP-TUV image and what is in it (section 5.1) | Only approved software may run | 3.4.8, 3.4.9 |
| 2 | **Allow rootless containers**: user namespaces must be enabled for the model user(s) | The RHEL 9 STIG sets `user.max_user_namespaces = 0` unless a documented need exists; rootless Podman needs them | 3.4.2, 3.4.6, 3.4.7 |
| 3 | **Allow outbound hosts**: for runs, the Copernicus Climate Data Store; for building on the VM, the build hosts (section 5.2) | Default-deny boundary | 3.13.1, 3.13.6 |
| 4 | **Crypto inside the container**: the only network traffic is the ERA5 download, made with RHEL 9's Python and OpenSSL, which follow the host's FIPS mode (section 5.5). The conda-forge libraries carry a second OpenSSL copy that the model never uses | FIPS-validated cryptography is required where crypto protects CUI | 3.13.11 |
| 5 | **The Copernicus API key**: where it is stored and who can read it | It is a credential for an external service | 3.1.1, 3.5.10 |
| 6 | **Whether the ASP decks and outputs are CUI**, and where they are stored | Storage, access and marking rules follow from this | 3.1.3, 3.8.x, 3.13.16 |
| 7 | **Image scanning and patching**: which scanner, how often the image is rebuilt | Vulnerability scanning and flaw remediation | 3.11.2, 3.11.3, 3.14.1 |
| 8 | **Long runs after logout** (`loginctl enable-linger`) | Runs take hours; idle sessions are locked or ended | 3.1.10, 3.1.11 |
| 9 | **Stack size limit**: an unlimited stack for the model user | The model needs it; some baselines lock limits | 3.4.2 |
| 10 | **No GPU access for the model's containers**: the model doesn't use the GPU, so its containers need no GPU device and the NVIDIA container toolkit isn't needed for it | Least functionality | 3.4.6, 3.4.7 |

On item 4: the ERA5 download carries only public data and the Copernicus
key. If the ISSO prefers the container to have no network at all, download
the ERA5 files outside the container and give them to the model (section
7.3); then run with `--network none`.

---

## 3. Machine size: CPU, memory, GPU

### 3.1 CPU (16 AMD vCPUs)

- The image runs **16 MPI ranks** by default, one per vCPU. Check what the
  vCPUs are: `lscpu` shows `Thread(s) per core`. If the hypervisor gives
  the VM 8 physical cores with 2 threads each, 16 ranks still work but each
  pair shares a core; try a build with `--build-arg NSUBDOMAINS_X=8` and
  compare the run time of a short case. Ask the virtualization administrator
  not to overcommit the VM's vCPUs: the model is CPU-bound.
- **AMD processors are fine.** The model is compiled with the Intel `ifx`
  compiler but without Intel-only instructions (no `-x` options), so it runs
  on any x86-64 processor. The Intel math run-time library picks its code
  path by processor, so results are bit-identical from run to run on the same
  machine, but the last digits can differ from a run on an Intel processor.
  Keep that in mind when comparing with results from another machine.
- Run **one case at a time**: two runs at once share the 16 vCPUs and each
  takes about twice as long.

### 3.2 Memory (64 GB) and disk

- 64 GB is ample. A 72-step test used about 1.1 GB in total (measured); the
  container guide plans conservatively for up to 0.7 GB per rank (11–12 GB
  for 16 ranks) on long runs.
- Disk: about 15 GB free to build the image (the Intel compiler alone is
  about 4 GB in the build stage), 0.75 GB for the finished image, and roughly
  1–2 GB of output per model hour. Rootless Podman keeps images in the
  user's home directory (`~/.local/share/containers`), so check
  home-directory quotas.

### 3.3 The NVIDIA L4 GPU is not used

SAM-ASP-TUV runs on CPUs only: it is Fortran parallelized with MPI, and has
no GPU code. The image therefore contains no CUDA, and its containers need
no GPU: don't install the NVIDIA container toolkit or pass `--device
nvidia.com/gpu=...` for this model (item 10). The GPU stays free for other
work on the VM.

Porting the model to GPUs would be a separate, large project (the Intel
compiler's GPU offload targets Intel GPUs; NVIDIA GPUs would need OpenACC or
CUDA Fortran with NVIDIA's compiler). The L4 would also be a weak fit: ASP's
chemistry is computed in double precision, and the L4's double-precision
throughput is a small fraction (1/64) of its single-precision throughput.

---

## 4. Administrator setup (once)

These steps need root and are done by the VM administrator.

**4.1 Check the baseline.**

```bash
cat /etc/os-release            # RHEL 9.x (or Rocky/Alma 9)
fips-mode-setup --check        # "FIPS mode is enabled."
getenforce                     # Enforcing
lscpu | grep -E '^CPU\(s\)|Thread|Core'
```

**4.2 Install Podman** (from the RHEL repositories, so it is covered by the
OS's patching):

```bash
dnf install -y container-tools        # podman, buildah, skopeo
```

**4.3 Enable rootless containers for the model user(s)** (after the ISSO
approves item 2). Replace `samuser` with the account name.

```bash
# user namespaces (the STIG may have set this to 0)
echo 'user.max_user_namespaces = 15000' > /etc/sysctl.d/90-samasp-userns.conf
sysctl --system

# subordinate UID/GID ranges; useradd adds them for local accounts,
# directory (IdM/AD) accounts need them added by hand
grep samuser /etc/subuid /etc/subgid || usermod --add-subuids 200000-265535 --add-subgids 200000-265535 samuser

# let the user's runs continue after logout (item 8)
loginctl enable-linger samuser
```

**4.4 Allow an unlimited stack** for the model user (item 9):

```bash
echo 'samuser  hard  stack  unlimited' > /etc/security/limits.d/90-samasp.conf
```

(The user then passes `--ulimit stack=-1` to Podman; a rootless container
can't raise a limit above the user's hard limit.)

**4.5 Proxy.** If the enclave uses a proxy, make sure the model user's
environment sets it (for example in `/etc/profile.d/proxy.sh`):

```bash
export https_proxy=http://proxy.example.org:3128  HTTPS_PROXY=$https_proxy
export no_proxy=localhost,127.0.0.1
```

Podman passes these variables into containers by default. If the proxy
inspects TLS, the organization's CA certificate must be in the host's trust
store (`/etc/pki/ca-trust/source/anchors/`, then `update-ca-trust`); the user
mounts it into the container (section 6).

**4.6 Storage for outputs.** Create the directory where model outputs will
live, on storage approved for CUI (encrypted at rest with FIPS-validated
encryption, 3.13.16; backed up per your policy), owned by the model user:

```bash
mkdir -p /data/samasp/samuser && chown samuser: /data/samasp/samuser && chmod 700 /data/samasp/samuser
```

Rootless Podman can't use network filesystems (NFS) for its own image
storage, but output directories on NFS work if the user can write to them.

---

## 5. Getting the image

### 5.1 What is in the image (for software approval)

The RHEL 9 image is built from `Dockerfile.ubi9`. Its parts:

| Component | Source |
|---|---|
| Base OS | Red Hat Universal Base Image 9.8 (`registry.access.redhat.com/ubi9/ubi:9.8`), or your organization's hardened UBI 9 image (build option `BASE`) |
| Python 3.9, NumPy, `requests`, OpenSSL | RHEL 9 packages from the UBI repositories (`cdn-ubi.redhat.com`) |
| Open MPI 4.1.6, HDF5 1.14.6, NetCDF-C 4.9.3 | conda-forge (`conda.anaconda.org`), pinned versions. They are not in the UBI repositories: Open MPI is in full RHEL, HDF5 and NetCDF only in EPEL |
| Intel Fortran/C compiler 2026.2.0 (run-time libraries only in the final image) | conda-forge, pinned; Intel's licence terms apply |
| NetCDF-Fortran 4.6.2 | Built from source with the Intel compiler; the source tarball is checked against a SHA-256 checksum |
| `cdsapi` 0.7.7, `netCDF4` 1.6.5, `cftime` 1.6.3 (and their dependencies) | PyPI (`pypi.org`, `files.pythonhosted.org`), pinned top-level versions |
| SAM-ASP-TUV, ASP (common_ASP), TUV | This repository and its submodule, compiled from source |
| Post-processing tools (`UTIL/`) | Compiled from source in this repository |

The image holds no CUI and no credentials: the Copernicus key and your ASP
decks are mounted at run time. It runs no services, opens no ports, and
contains no CUDA or GPU software.

To list what is installed (for the configuration records or a scan review):

```bash
podman run --rm --entrypoint rpm sam-asp-tuv -qa | sort                  # RHEL packages
podman run --rm --entrypoint ls sam-asp-tuv /opt/deps/conda-meta         # conda-forge packages
podman run --rm --entrypoint python3 sam-asp-tuv -c \
  "import importlib.metadata as m; print(*sorted(f'{d.name} {d.version}' for d in m.distributions()), sep='\n')"
```

### 5.2 Option A: build on the VM

This needs these hosts (item 3):

| Host | For |
|---|---|
| `github.com` | the code (or copy the repository in) |
| `registry.access.redhat.com`, and the download host it redirects image layers to (`cdn01.quay.io` when this was tested) | the UBI 9 base image (not needed if `BASE` points to an internal registry) |
| `cdn-ubi.redhat.com` | UBI 9 packages |
| `conda.anaconda.org` | Intel compiler, Open MPI, HDF5, NetCDF-C |
| `archive.ubuntu.com` | the NetCDF-Fortran source tarball only. To use the upstream release or an internal mirror instead, pass `--build-arg NETCDF_FORTRAN_URL=...` and the matching `--build-arg NETCDF_FORTRAN_SHA256=...` |
| `pypi.org`, `files.pythonhosted.org` | `cdsapi`, `netCDF4` |

If a download fails, the proxy logs show which host was refused.

As the model user:

```bash
git clone https://github.com/matt-alvarado/SAM-ASP-TUV.git
cd SAM-ASP-TUV
git checkout asp_v2.2.2
git submodule update --init

podman build -f Dockerfile.ubi9 -t sam-asp-tuv:$(git rev-parse --short HEAD) \
  --secret id=ca,src=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem \
  .
podman tag sam-asp-tuv:$(git rev-parse --short HEAD) sam-asp-tuv:latest
```

- To build on your organization's hardened UBI 9 image, add `--build-arg
  BASE=<registry>/<hardened-ubi9-image>:<tag>`. It must still be able to
  install the RHEL packages listed above with `dnf` (from the UBI
  repositories or your internal mirror).
- The `--secret` line gives the build the host's trust store, so downloads
  work through a TLS-inspecting proxy. The certificates are used only during
  the build and are not left in the image.
- Podman passes the host's proxy variables to the build.
- The build takes 15–30 minutes and needs about 15 GB of free disk (3.2).

### 5.3 Option B: build elsewhere and transfer the image

If the ISSO prefers not to open the build hosts from the enclave, build the
image with the same `podman build -f Dockerfile.ubi9 ...` command on a
machine outside the enclave (any x86-64 Linux machine with Podman, or Docker
with BuildKit), then:

```bash
# outside the enclave
podman save -o sam-asp-tuv-<version>.tar sam-asp-tuv:<version>
sha256sum sam-asp-tuv-<version>.tar > sam-asp-tuv-<version>.tar.sha256
```

Bring both files in through your approved transfer process (3.8.7, 3.8.9),
then on the VM:

```bash
sha256sum -c sam-asp-tuv-<version>.tar.sha256
podman load -i sam-asp-tuv-<version>.tar
podman tag sam-asp-tuv:<version> sam-asp-tuv:latest
```

Images saved by Docker load into Podman the same way.

### 5.4 Scan the image

Before first use and after every rebuild, scan the image with the
vulnerability scanner your organization uses (3.11.2), and record the result
with the image ID (`podman image inspect --format '{{.Id}}' sam-asp-tuv`).

### 5.5 Check FIPS mode inside the container

On a FIPS-mode RHEL 9 host, RHEL 9 containers run their OpenSSL in FIPS mode
too: they read the host kernel's FIPS setting. Check it once on the VM:

```bash
podman run --rm --entrypoint cat sam-asp-tuv /proc/sys/crypto/fips_enabled   # 1
podman run --rm --entrypoint openssl sam-asp-tuv list -providers             # lists a "fips" provider
podman run --rm --entrypoint python3 sam-asp-tuv -c "import ssl; print(ssl.OPENSSL_VERSION)"
```

The Python used for the ERA5 download is RHEL's, linked to RHEL's OpenSSL
(`/usr/lib64/libssl.so.3`). The conda-forge NetCDF-C library brings its own
OpenSSL copy (`/opt/deps/lib`) for its remote-data features, and the
`netCDF4` Python package bundles one too; the model reads only local files,
so neither is used. Note both in the software record (item 4).

---

## 6. User setup (once per user)

As the model user (no root needed):

**6.1 Check Podman works rootless:**

```bash
podman info --format '{{.Host.Security.Rootless}} {{.Host.IDMappings.UIDMap}}'
ulimit -Hs          # should print "unlimited"
```

**6.2 Set up ERA5 access.** Create the Copernicus account, accept the two
ERA5 licences and get the API key as in
[2.1 of the container guide](Container_User_Guide.md#21-setting-up-era5-access-once).
Store the key in a file only you can read, not in an environment variable
(environment variables show in `podman inspect` and process listings):

```bash
umask 077
cat > ~/.cdsapirc <<'EOF'
url: https://cds.climate.copernicus.eu/api
key: <your personal access token>
EOF
chmod 600 ~/.cdsapirc
```

**6.3 Organization CA certificate** (only if the proxy inspects TLS). Copy
the host's trust store to your own directory, so it can be mounted over the
container's trust store without relabelling a system file:

```bash
mkdir -p ~/samasp
cp /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem ~/samasp/ca-bundle.pem
```

**6.4 Directories.** Keep your ASP decks and outputs on the approved storage
from 4.6, for example:

```bash
D=/data/samasp/$USER
mkdir -p $D/decks $D/runs
cp -r /path/to/your/decks/* $D/decks/      # the seven .in files
```

---

## 7. Running a case

### 7.1 One case

```bash
D=/data/samasp/$USER
mkdir -p $D/runs/williams

podman run --rm --ulimit stack=-1 --shm-size=1g \
  -v $D/decks:/decks:ro,Z \
  -v $D/runs/williams:/output:Z \
  -v ~/.cdsapirc:/root/.cdsapirc:ro,Z \
  -v ~/samasp/ca-bundle.pem:/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:ro,Z \
  -e ASP_DECK_DIR=/decks -e CASE=williams \
  -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 \
  -e RUN_HOURS=4.5 \
  sam-asp-tuv
```

- `:Z` gives the container (and only that container) SELinux access to the
  directory or file. Use it only on your own files and directories, never on
  system paths or your whole home directory.
- Leave out the `ca-bundle.pem` line if the proxy doesn't inspect TLS.
- Root inside a rootless container is your own account outside it, so the
  outputs belong to you and the container has no more rights than you do.
- No `--device` or GPU option: the model doesn't use the GPU (3.3).
- All settings (`LAT`, `RUN_HOURS`, `FIRE_START_STEP`, ...) are in
  [section 5 of the container guide](Container_User_Guide.md#5-settings-reference).

### 7.2 Long runs: run it as a systemd user service

Runs take hours, and the enclave's session policy will lock or end your
login session. Start the run as a systemd user service instead; with
lingering enabled (4.3) it keeps going after you log out, and its output
goes to the system journal (3.3.1):

```bash
systemd-run --user --unit=samasp-williams --collect \
  podman run --rm --ulimit stack=-1 --shm-size=1g \
    -v $D/decks:/decks:ro,Z -v $D/runs/williams:/output:Z \
    -v $HOME/.cdsapirc:/root/.cdsapirc:ro,Z \
    -v $HOME/samasp/ca-bundle.pem:/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:ro,Z \
    -e ASP_DECK_DIR=/decks -e CASE=williams \
    -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 \
    -e RUN_HOURS=4.5 \
    sam-asp-tuv
```

Then:

```bash
systemctl --user status samasp-williams        # running?
journalctl --user -u samasp-williams -f        # follow the container's messages
tail -f $D/runs/williams/run.out               # follow the model log
systemctl --user stop samasp-williams          # stop the run
```

### 7.3 Without network access from the container

If the ISSO decides the container should have no network (item 4), download
the two ERA5 files elsewhere (see
[7.3 of the container guide](Container_User_Guide.md#73-using-era5-files-you-already-have)),
put them on the VM, and run with `--network none`:

```bash
podman run --rm --network none --ulimit stack=-1 --shm-size=1g \
  -v $D/decks:/decks:ro,Z -v $D/era5:/era5:ro,Z -v $D/runs/williams:/output:Z \
  -e ASP_DECK_DIR=/decks -e CASE=williams \
  -e LAT=39.70 -e LON=-122.58 -e DATE=2009-11-17 -e START_TIME=20:00 -e RUN_HOURS=4.5 \
  -e ERA5_PL_FILE=/era5/pl.nc -e ERA5_SL_FILE=/era5/sl.nc \
  sam-asp-tuv
```

No key or CA certificate is needed then.

### 7.4 Converting outputs

The conversion tools run the same way, for example:

```bash
podman run --rm --network none -v $D/runs/williams/OUT_STAT:/data:Z -w /data \
  --entrypoint /opt/sam-asp-tuv/UTIL/stat2nc sam-asp-tuv williams_williams.stat
```

---

## 8. Handling inputs and outputs

- **ASP decks.** If your decks hold CUI (for example, emission data with
  distribution limits), keep them on the approved storage only, and never
  build them into an image: an image with CUI in it is itself CUI and can't
  be built or stored outside the enclave.
- **Outputs.** Treat the run directory under the rules for your project's
  data: mark it as your CUI program requires, keep it on the approved
  storage, and move it out of the enclave only through approved processes.
  `case/met_source.txt` and `ERA5_cache/` hold only public ERA5 data.
- **Access.** Keep run directories readable only by the people who need them
  (`chmod 700`, or a project group), 3.1.1–3.1.3.
- **The Copernicus key** never goes into the run directory, logs, scripts
  or the image.

---

## 9. Records, updates and removal

**For each run**, keep (this also makes runs repeatable; see
[13 of the container guide](Container_User_Guide.md#13-reproducibility)):
the image ID and tag, the exact `podman run` command, the ASP decks, and
the run directory's `case/`, `ERA5_cache/` and `run.out`.

```bash
podman image inspect --format '{{.Id}} {{.RepoTags}}' sam-asp-tuv
```

**Updates.** The image's RHEL packages get security fixes over time (a
rebuild installs the current ones). Rebuild the image on the schedule your
flaw-remediation policy sets (for example monthly, or when a scan finds a
fixed vulnerability; 3.14.1), scan it again, and keep the previous image
until results from the new one are checked. The compiler, MPI, NetCDF and
Python packages are pinned in `Dockerfile.ubi9`; to take a fix in one of
them, change its version there. Rebuilding with the same code and pinned
versions should give the same model results; note the image ID with every
run.

**Removal.** To remove the model from the VM:

```bash
podman rmi --all            # as each model user
podman system reset         # also removes the user's container storage
```

and delete the run directories as your media-sanitization procedure
requires (3.8.3). The administrator removes `container-tools` and the
settings from section 4 if nothing else needs them.

---

## 10. Troubleshooting

| Problem | Cause and fix |
|---|---|
| `cannot find newuidmap` or `there might not be enough IDs available in the namespace` | The user has no subordinate UID/GID range, or `newuidmap` isn't installed: section 4.3. After adding ranges, run `podman system migrate`. |
| `user namespaces are not enabled` / `operation not permitted` at start | `user.max_user_namespaces` is 0 (STIG setting): section 4.3, after ISSO approval. |
| `Permission denied` reading decks or writing outputs | SELinux: add `:Z` to the `-v` option (section 7.1). Check `ausearch -m avc -ts recent` (administrator). Also check the directory's ordinary permissions. |
| `WARNING: could not lift the stack limit` | The user's hard stack limit is not unlimited: section 4.4, then log in again. |
| The run stops when you log out | Lingering isn't enabled, or the run wasn't started with `systemd-run`: sections 4.3 and 7.2. |
| ERA5 download: `SSL: CERTIFICATE_VERIFY_FAILED` | The proxy inspects TLS: mount the CA bundle (section 6.3, 7.1). |
| ERA5 download: connection refused or timed out | The proxy variables aren't set, or the Copernicus hosts aren't allowed. Copernicus serves finished downloads from a separate host; the proxy logs show which host was refused. |
| `podman build`: `short-name resolution` prompt or error | A `BASE` given without a registry name; give the full name (`registry.example.org/...`). |
| `podman build`: `sha256sum: WARNING: computed checksum did NOT match` | `NETCDF_FORTRAN_URL` was changed without the matching `NETCDF_FORTRAN_SHA256`, or the download was altered. |
| `run.out`: `plm_rsh_agent ... ssh : rsh` could not be found | Open MPI looked for `ssh`. The RHEL 9 image sets `OMPI_MCA_plm=isolated` (one machine, no `ssh`); don't override it with `-e OMPI_MCA_plm=...`. |
| `podman build`: download errors | Build hosts not allowed (section 5.2), or the CA bundle wasn't passed with `--secret`. |
| Image storage fills the home directory | Ask the administrator for space, or set `graphroot` in `~/.config/containers/storage.conf` to a larger local (not NFS) directory you own. |

Model problems (missing deck files, ERA5 coverage, fire timing) are covered
in [section 12 of the container guide](Container_User_Guide.md#12-troubleshooting).

---

## 11. What was tested

In this project's cloud workspace (Ubuntu 24.04 host, Podman 4.9):

- **The RHEL 9 image** (`Dockerfile.ubi9` on UBI 9.8): built with `podman
  build` through a TLS-inspecting proxy, with the CA bundle as a build
  secret. A 72-step fire case with mounted decks and ERA5 files, run with
  `--network none`, finished, and its statistics and 2-D output files were
  **byte-identical** to the same run with the Ubuntu-based image.
- **Rootless Podman** as an ordinary user, with `:Z` mounts: the run
  finished and the outputs belonged to the user (see the note below on which
  image).
- **Python** uses RHEL's OpenSSL 3.5.8 inside the image.
- **Memory**: the 72-step test peaked at about 1.1 GB for 16 ranks.

Not tested: a RHEL 9 host, FIPS mode (section 5.5), the STIG settings,
SELinux enforcing, an organization's hardened UBI 9 base, `systemd-run` with
lingering, an enterprise proxy at run time, real ERA5 downloads, and the
AMD/L4 VM itself. The steps follow Red Hat's and Podman's documentation; try
a short case first (for example `RUN_HOURS=0.2 FIRE_START_STEP=0`) before
planning long runs.
