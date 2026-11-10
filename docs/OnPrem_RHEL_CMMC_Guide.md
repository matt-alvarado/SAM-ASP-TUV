# SAM-ASP-TUV on an on-prem RHEL VM in a CMMC Level 2 enclave

This guide covers installing and running the SAM-ASP-TUV container on an
on-premises Linux virtual machine (VM) inside a CMMC Level 2 enclave, that
is, a system that handles Controlled Unclassified Information (CUI) and must
meet the 110 requirements of NIST SP 800-171 rev 2.

It adds what is specific to that setting to the
[container user's guide](Container_User_Guide.md), which explains the model
settings, ASP decks, ERA5, fire timing and outputs. Read that guide for
everything about the model itself.

> **This is not a compliance determination.** CMMC sets requirements, not
> products. What is allowed on your VM is decided by your organization's
> System Security Plan (SSP) and its security officer (ISSO/ISSM). The
> control numbers below (e.g. 3.4.8) point to the NIST SP 800-171
> requirements each step touches, so you and the ISSO can find them. Where
> your SSP says otherwise, follow your SSP.

Contents:
1. [Assumptions](#1-assumptions)
2. [What to agree with your ISSO first](#2-what-to-agree-with-your-isso-first)
3. [Machine size](#3-machine-size)
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
| VM | Inside the CUI enclave; 16 cores, 64 GB RAM (4 GB per core) |
| OS | Red Hat Enterprise Linux 9 (or Rocky/Alma 9), hardened (DISA STIG or CIS), **FIPS mode on**, SELinux enforcing |
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
| 1 | **Approve the software**: Podman (`container-tools`) on the VM, and the SAM-ASP-TUV image and what is in it (section 5.1) | Only approved software may run | 3.4.8, 3.4.9 |
| 2 | **Allow rootless containers**: user namespaces must be enabled for the model user(s) | The RHEL 9 STIG sets `user.max_user_namespaces = 0` unless a documented need exists; rootless Podman needs them | 3.4.2, 3.4.6, 3.4.7 |
| 3 | **Allow outbound hosts**: for runs, the Copernicus Climate Data Store; for building on the VM, the build hosts (section 5.2) | Default-deny boundary | 3.13.1, 3.13.6 |
| 4 | **Crypto inside the container**: the image is Ubuntu-based, so its TLS (used only to download public ERA5 data from Copernicus) does not use FIPS-validated modules, even on a FIPS-mode host | FIPS-validated cryptography is required where crypto protects CUI | 3.13.11 |
| 5 | **The Copernicus API key**: where it is stored and who can read it | It is a credential for an external service | 3.1.1, 3.5.10 |
| 6 | **Whether the ASP decks and outputs are CUI**, and where they are stored | Storage, access and marking rules follow from this | 3.1.3, 3.8.x, 3.13.16 |
| 7 | **Image scanning and patching**: which scanner, how often the image is rebuilt | Vulnerability scanning and flaw remediation | 3.11.2, 3.11.3, 3.14.1 |
| 8 | **Long runs after logout** (`loginctl enable-linger`) | Runs take hours; idle sessions are locked or ended | 3.1.10, 3.1.11 |
| 9 | **Stack size limit**: an unlimited stack for the model user | The model needs it; some baselines lock limits | 3.4.2 |

On item 4: the only network traffic from the container is the ERA5 download
from Copernicus (public data; the request carries the Copernicus key). If
your ISSO requires FIPS-validated cryptography inside containers anyway, the
image can be rebuilt on Red Hat's Universal Base Image (UBI 9); ask the
model maintainers. Alternatively, download the ERA5 files outside the
container and give them to the model ([7.3 of the container guide](Container_User_Guide.md#73-using-era5-files-you-already-have)),
so the container needs no network at all (`--network none`).

---

## 3. Machine size

With 16 cores and 64 GB:

- **CPU.** The image runs 16 MPI ranks by default, one per core, which fits
  this VM. Check what the VM's cores are: `lscpu` shows `Thread(s) per
  core`. Ask the virtualization administrator not to overcommit the VM's
  vCPUs; the model is CPU-bound, and shared or hyper-threaded vCPUs make it
  slower. Run **one case at a time**: two runs at once share the 16 cores
  and each takes about twice as long.
- **Memory.** 64 GB is ample. A 72-step test used about 1.1 GB in total
  (measured); the container guide plans conservatively for up to 0.7 GB per
  rank (11–12 GB for 16 ranks) on long runs.
- **Disk.** About 10 GB to build the image (if you build on the VM), 1 GB
  for the image itself, and roughly 1–2 GB of output per model hour. Rootless
  Podman keeps images in the user's home directory
  (`~/.local/share/containers`), so check home-directory quotas.

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

| Component | Source |
|---|---|
| Base OS | Ubuntu 24.04 (`docker.io/library/ubuntu:24.04`) |
| MPI, NetCDF-C, Python 3, NumPy, netCDF4 | Ubuntu 24.04 packages (`archive.ubuntu.com`, `security.ubuntu.com`) |
| Intel Fortran/C compiler run-time libraries (`ifx`, `icx`) | conda-forge (`conda.anaconda.org`); Intel's licence terms apply |
| NetCDF-Fortran | Built from Ubuntu's source package with the Intel compiler |
| `cdsapi` (Copernicus client) | PyPI (`pypi.org`, `files.pythonhosted.org`) |
| SAM-ASP-TUV, ASP (common_ASP), TUV | This repository and its submodule, compiled from source |
| Post-processing tools (`UTIL/`) | Compiled from source in this repository |

The image holds no CUI and no credentials: the Copernicus key and your ASP
decks are mounted at run time. It runs no services and opens no ports.

To list the installed packages (for the configuration records or a scan
review):

```bash
podman run --rm --entrypoint dpkg-query sam-asp-tuv -W -f '${Package} ${Version}\n'
podman run --rm --entrypoint python3 sam-asp-tuv -c \
  "import importlib.metadata as m; print(*sorted(f'{d.name} {d.version}' for d in m.distributions()), sep='\n')"
```

### 5.2 Option A: build on the VM

This needs the build hosts in item 3: `github.com` (the code; or copy the
repository in), the Docker Hub registry (`registry-1.docker.io`,
`auth.docker.io`, `production.cloudflare.docker.com`),
`archive.ubuntu.com`, `security.ubuntu.com`, `conda.anaconda.org` and
`pypi.org`/`files.pythonhosted.org`. If a download fails, the proxy logs
show which host it was refused.

As the model user:

```bash
git clone https://github.com/matt-alvarado/SAM-ASP-TUV.git
cd SAM-ASP-TUV
git checkout asp_v2.2.2
git submodule update --init

podman build -t sam-asp-tuv:$(git rev-parse --short HEAD) \
  --secret id=ca,src=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem \
  .
podman tag sam-asp-tuv:$(git rev-parse --short HEAD) sam-asp-tuv:latest
```

The `--secret` line gives the build the host's trust store, so downloads
work through a TLS-inspecting proxy. The certificates are used only during
the build and are not left in the image. Podman passes the host's proxy
variables to the build. The build takes 15–30 minutes.

### 5.3 Option B: build elsewhere and transfer the image

If the ISSO prefers not to open the build hosts from the enclave, build the
image on a machine outside the enclave (any Linux machine with Podman or
Docker; see [section 3 of the container guide](Container_User_Guide.md#3-building-the-image)),
then:

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
the host's trust store to your own directory, so it can be mounted into the
container without relabelling a system file:

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
  -v ~/samasp/ca-bundle.pem:/etc/ssl/certs/ca-certificates.crt:ro,Z \
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
    -v $HOME/samasp/ca-bundle.pem:/etc/ssl/certs/ca-certificates.crt:ro,Z \
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

**Updates.** The image's Ubuntu packages get security fixes over time.
Rebuild the image on the schedule your flaw-remediation policy sets (for
example monthly, or when a scan finds a fixed vulnerability; 3.14.1), scan
it again, and keep the previous image until results from the new one are
checked. Rebuilding with the same code gives the same model results unless
the compiler version changed, so note the image ID with every run.

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
| `podman build`: `short-name resolution` prompt or error | An old copy of the Dockerfile names the base image `ubuntu:24.04`; use `--build-arg UBUNTU=docker.io/library/ubuntu:24.04`. |
| `podman build`: download errors | Build hosts not allowed (section 5.2), or the CA bundle wasn't passed with `--secret`. |
| Image storage fills the home directory | Ask the administrator for space, or set `graphroot` in `~/.config/containers/storage.conf` to a larger local (not NFS) directory you own. |

Model problems (missing deck files, ERA5 coverage, fire timing) are covered
in [section 12 of the container guide](Container_User_Guide.md#12-troubleshooting).

---

## 11. What was tested

- **Rootless Podman** (Podman 4.9 on Ubuntu 24.04, as an ordinary user):
  loading a Docker-built image, and a 72-step fire case with mounted decks
  and ERA5 files, `:Z` mounts and `--network none`. The run finished, the
  outputs belonged to the user, and the statistics and 2-D output files were
  byte-identical to the same run under Docker.
- **`podman build`** of the Dockerfile (Podman 4.9, run as root here, through
  a TLS-inspecting proxy with the CA bundle passed as a build secret): the
  build succeeded, no CA certificate was left in the image, and the 72-step
  case gave byte-identical results to the Docker-built image. A rootless
  `podman build` was not tried.
- **Memory**: the 72-step test peaked at about 1.1 GB for 16 ranks.

Not tested: RHEL 9 itself, FIPS mode, the STIG settings, SELinux enforcing,
`systemd-run` with lingering, a TLS-inspecting enterprise proxy at run time,
and real ERA5 downloads. The steps follow Red Hat's and Podman's
documentation; try a short case first (for example `RUN_HOURS=0.2
FIRE_START_STEP=0`) before planning long runs.
