# SAM-ASP-TUV on an Apple silicon Mac in a CMMC Level 2 enclave

This guide covers running SAM-ASP-TUV with **rootless Podman** on an **Apple
silicon (ARM64) Mac** that handles Controlled Unclassified Information (CUI)
and so must meet CMMC Level 2 (NIST SP 800-171 rev 2).

It builds on two other guides:

- [Container user's guide](Container_User_Guide.md): model settings, ASP
  decks, ERA5, outputs.
- [On-prem RHEL guide](OnPrem_RHEL_CMMC_Guide.md): the CMMC items, rootless
  Podman set-up, data handling and records for a RHEL 9 host. Most of it
  applies unchanged inside the Mac's RHEL 9 virtual machine; this guide gives
  what is different on a Mac.

> **This is not a compliance determination.** What is allowed on the Mac is
> decided by your System Security Plan (SSP) and your security officer
> (ISSO). Requirement numbers (e.g. 3.13.11) point to NIST SP 800-171.

Contents:
1. [How it fits together](#1-how-it-fits-together)
2. [What to agree with your ISSO](#2-what-to-agree-with-your-isso)
3. [The Mac](#3-the-mac)
4. [The RHEL 9 virtual machine](#4-the-rhel-9-virtual-machine)
5. [Building the ARM64 image](#5-building-the-arm64-image)
6. [Running cases](#6-running-cases)
7. [The GNU (gfortran) build](#7-the-gnu-gfortran-build)
8. [Alternative: Podman Desktop / podman machine](#8-alternative-podman-desktop--podman-machine)
9. [Troubleshooting](#9-troubleshooting)
10. [What was tested](#10-what-was-tested)

---

## 1. How it fits together

Containers are a Linux feature, so on macOS they always run inside a Linux
virtual machine (VM):

```
Mac (macOS, managed, FileVault)
 └─ hypervisor (Apple Virtualization framework: UTM, VMware Fusion or Parallels)
     └─ RHEL 9 for ARM (aarch64) VM: FIPS mode, STIG-hardened, SELinux enforcing
         └─ rootless Podman (an ordinary user, no root)
             └─ SAM-ASP-TUV container: UBI 9 (RHEL 9) ARM64, built with gfortran
```

- **The image** is built from `Dockerfile.ubi9-gnu`. The model normally uses
  Intel's `ifx` compiler, which produces only x86-64 code, so the ARM64 image
  is built with the GNU compilers (gfortran/gcc). It runs natively at full
  speed; see [section 7](#7-the-gnu-gfortran-build) for how its results
  compare.
- **The VM** is RHEL 9, as on the enclave's servers, so the
  [on-prem guide](OnPrem_RHEL_CMMC_Guide.md) (administrator set-up, rootless
  Podman, FIPS check, data handling, records) applies inside it.

## 2. What to agree with your ISSO

In addition to the items in [section 2 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#2-what-to-agree-with-your-isso-first):

| # | Item | Why | Requirement |
|---|---|---|---|
| M1 | **The Mac is in the CUI boundary**: enrolled in your MDM, FileVault on, supported macOS kept patched, screen lock, firewall, a non-admin daily account | Endpoint protection for CUI; macOS's encryption (corecrypto) has FIPS 140 validations | 3.1.10, 3.4.1-3.4.2, 3.13.16, 3.14.1 |
| M2 | **The hypervisor** (UTM, VMware Fusion or Parallels) is approved software | Only approved software may run | 3.4.8, 3.4.9 |
| M3 | **The guest OS is RHEL 9 for ARM in FIPS mode**, and Red Hat's FIPS validation covers the ARM64 (aarch64) platform you use: ask Red Hat or check the CMVP certificates | FIPS-validated cryptography where crypto protects CUI | 3.13.11 |
| M4 | **VM isolation settings**: shared folders only for the approved data directory, clipboard and drag-and-drop sharing off, USB pass-through off, NAT networking | Limit information flow between macOS and the VM | 3.1.3, 3.13.1 |
| M5 | **VM disk and snapshots** hold CUI: they stay on the FileVault volume; backups of them (e.g. Time Machine) are encrypted and approved; snapshots are deleted per your media policy | Media protection | 3.8.1, 3.8.3, 3.8.9 |
| M6 | **Results from the gfortran build** are accepted for your use (they are close to, but not bit-identical with, the Intel builds; [section 7](#7-the-gnu-gfortran-build)) | Configuration records | 3.4.1 |

## 3. The Mac

For an M-series **Pro or Max with 32 GB or more**:

| | Recommendation |
|---|---|
| Performance cores | Pro chips have 5-10, Max chips 10-12 (`sysctl hw.perflevel0.physicalcpu` shows the number). Give the model at most that many MPI ranks: macOS may run extra VM threads on the slower efficiency cores. |
| MPI ranks | Build with `NSUBDOMAINS_X=8` (must divide 240; 10 and 12 also work on 10- or 12-core Max chips). |
| VM | CPUs = number of performance cores (e.g. 10), memory 16 GB (24 GB with 64 GB in the Mac), disk 60 GB (the build needs about 10 GB; outputs about 1-2 GB per model hour). |
| Power | Keep the Mac on power and stop it sleeping during runs (`caffeinate -dimsu` in a macOS terminal, or the Energy settings); runs take hours. |

A 72-step test with 16 ranks used about 1.1 GB of memory, so memory is not the
limit; cores are.

## 4. The RHEL 9 virtual machine

1. **Hypervisor.** Install the approved hypervisor (item M2). All three listed
   run ARM64 Linux guests on Apple silicon.
2. **RHEL 9 for ARM.** Install "Red Hat Enterprise Linux 9 for ARM 64"
   (aarch64) from your organization's subscription, with the hardened profile
   your SSP requires (DISA STIG or CIS). Enable FIPS mode **during
   installation** (installer boot option `fips=1`), which is the supported
   way, then check after the first boot:

   ```bash
   fips-mode-setup --check          # FIPS mode is enabled.
   uname -m                         # aarch64
   getenforce                       # Enforcing
   ```
3. **VM settings** (item M4): the CPU, memory and disk from
   [section 3](#3-the-mac); NAT networking; clipboard, drag-and-drop and USB
   sharing off; at most one shared folder, for approved data (or none: copy
   files with `scp` to the VM instead).
4. **Inside the VM**, follow the on-prem guide:
   [4 Administrator setup](OnPrem_RHEL_CMMC_Guide.md#4-administrator-setup-once)
   (Podman, user namespaces, subordinate IDs, lingering, stack limit, proxy,
   output storage) and
   [6 User setup](OnPrem_RHEL_CMMC_Guide.md#6-user-setup-once-per-user)
   (Copernicus key file, CA certificate, directories).

## 5. Building the ARM64 image

Inside the VM, as the model user:

```bash
git clone https://github.com/matt-alvarado/SAM-ASP-TUV.git
cd SAM-ASP-TUV
git checkout asp_v2.2.2
git submodule update --init

podman build -f Dockerfile.ubi9-gnu --build-arg NSUBDOMAINS_X=8 \
  -t sam-asp-tuv:gnu-8r-$(git rev-parse --short HEAD) \
  --secret id=ca,src=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem \
  .
podman tag sam-asp-tuv:gnu-8r-$(git rev-parse --short HEAD) sam-asp-tuv
```

- It builds natively for ARM64 (Podman builds for the VM's own architecture).
- Hosts it needs (item 3 of the on-prem guide): `github.com`,
  `registry.access.redhat.com` (and the host it redirects image layers to,
  `cdn01.quay.io` when this was tested), `cdn-ubi.redhat.com`,
  `conda.anaconda.org`, `pypi.org` and `files.pythonhosted.org`.
- `--build-arg BASE=<registry>/<hardened-ubi9>:<tag>` uses your organization's
  hardened UBI 9 image, if it has an ARM64 variant.
- For a case with other decks, add `--build-arg NTRACERS=<n>` (for example
  1486 for `cases/williams_2009`); the entry point reports the value if it is
  wrong, before doing anything.
- Instead of building in the VM, the image can be built on another ARM64
  machine and transferred (`podman save` / `podman load`, as in
  [5.3 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#53-option-b-build-elsewhere-and-transfer-the-image)).
  An image built on an x86-64 machine with `--platform linux/arm64` (under
  emulation) also works but takes much longer to build.

Then scan it and check FIPS mode inside the container as in
[5.4-5.5 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#54-scan-the-image).

## 6. Running cases

Exactly as in [section 7 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#7-running-a-case),
inside the VM; for a prepared case:

```bash
D=/data/samasp/$USER
mkdir -p $D/runs/wf2019
systemd-run --user --unit=samasp-wf2019 --collect \
  podman run --rm --ulimit stack=-1 --shm-size=1g \
    --env-file $HOME/SAM-ASP-TUV/cases/williams_flats_2019/case.env \
    -v $D/runs/wf2019:/output:Z \
    -v $HOME/.cdsapirc:/root/.cdsapirc:ro,Z \
    -v $HOME/samasp/ca-bundle.pem:/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:ro,Z \
    sam-asp-tuv
```

Keep the Mac awake and on power for the length of the run (section 3).

## 7. The GNU (gfortran) build

The ARM64 image is the same model, decks and settings, compiled with gfortran
13.4 and gcc instead of Intel ifx/icx (`build_gnu.sh`, which is
`build_ifx.sh` with `SAMASP_FC=gfortran`; the Makefiles of SAM, ASP and TUV
and `UTIL/build_util.sh` have GNU settings).

- **Floating point.** Fused multiply-add contraction is turned off
  (`-ffp-contract=off`), so the ARM64 and x86-64 gfortran builds compute the
  same arithmetic; ASP uses the GNU settings of the ASP repository's own
  gfortran build (`-frounding-math -fsignaling-nans`, 8-byte reals).
- **Agreement with the Intel build** (72-step test case, x86-64): plume means
  of CO, O3, NO, NO2, PAN, HONO, OH, total and secondary organic aerosol and
  the photolysis rates agree within about 0.5%, the particle number of the
  main size bins within 0.1%. Point-by-point values differ more, because the
  large-eddy turbulence amplifies tiny arithmetic differences; that is
  expected of any change of compiler or machine.
- **Reproducibility.** The same image on the same machine and rank count gives
  the same results run after run. Results are **not** bit-identical with the
  Intel builds, or between a Mac and an x86-64 machine; compare runs made with
  the same build (item M6).

## 8. Alternative: Podman Desktop / podman machine

Podman Desktop (or `podman machine init`) is the usual way to run Podman on a
Mac. It creates its own VM running **Fedora CoreOS**, with rootless containers
by default. It is simpler, but **Fedora's cryptographic modules are not FIPS
140 validated**, so it does not meet item M3 when the VM or container uses
cryptography to protect CUI. Use it only if your ISSO accepts it, for example:

- for development and tests with non-CUI data; or
- with the container offline (`--network none`, ERA5 files supplied, see
  [7.3 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#73-without-network-access-from-the-container))
  and the VM disk on the FileVault volume.

```bash
podman machine init --cpus 10 --memory 16384 --disk-size 60   # rootless by default
podman machine start
podman build -f Dockerfile.ubi9-gnu --build-arg NSUBDOMAINS_X=8 -t sam-asp-tuv .
```

Files under your macOS home directory are shared into this VM by default;
restrict that to the data directory (`--volume` on `podman machine init`)
per item M4.

## 9. Troubleshooting

| Problem | Cause and fix |
|---|---|
| `exec format error` | An x86-64 image on the ARM64 VM: build with `Dockerfile.ubi9-gnu` in the VM, or load an ARM64 image. |
| `need N SAM tracers, but this model was built with ntracers = M` | Rebuild with `--build-arg NTRACERS=N`. Nothing was run. |
| `NSUBDOMAINS_X=... must divide nx_gl=240` | Use 4, 6, 8, 10, 12, 16, 20, ... |
| Runs much slower than expected | More ranks than performance cores, the Mac on battery or throttling, or the VM given fewer CPUs than ranks. |
| The run stops when the Mac sleeps | Use `caffeinate` (section 3) and `systemd-run` (section 6). |

Other problems: [section 10 of the on-prem guide](OnPrem_RHEL_CMMC_Guide.md#10-troubleshooting)
and [section 12 of the container guide](Container_User_Guide.md#12-troubleshooting).

## 10. What was tested

In this project's cloud workspace (x86-64 Linux, no Mac available):

- **gfortran build, x86-64, native**: model, ASP, TUV and all 23 `UTIL` tools
  build; the 72-step test case runs and agrees with the Intel build as in
  section 7. This comparison also found a bug in `write_fields3D.f90` (five
  gas-phase organic output fields read past their arrays), fixed for all builds.
- **`Dockerfile.ubi9-gnu`, x86-64**: builds and runs the test case.
- **`Dockerfile.ubi9-gnu`, ARM64**: see the note below on the build and test
  run under QEMU emulation.

Not tested: a Mac, the hypervisors, RHEL 9 for ARM, FIPS mode, Podman machine,
and run times on Apple silicon (emulation times say nothing about native
speed). Run a short case first (for example `--hours 0.2 --fire-start-step 0`)
and time it before planning long runs.
