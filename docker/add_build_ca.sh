#!/bin/sh
# Trust the extra CA certificate given to a container build as the secret
# "ca" (docker/podman build --secret id=ca,src=ca.crt), for building behind a
# TLS-intercepting proxy; does nothing without the secret. With --remove,
# drop it again (the runtime image must not keep it).
# Works on Debian/Ubuntu (update-ca-certificates; the ca-certificates package
# must be installed) and RHEL/UBI (update-ca-trust).
set -e
if [ -d /etc/pki/ca-trust/source/anchors ]; then
    dest=/etc/pki/ca-trust/source/anchors/build-ca.pem
    refresh() { update-ca-trust; }
else
    dest=/usr/local/share/ca-certificates/build-ca.crt
    refresh() { update-ca-certificates "$@" > /dev/null; }
fi
if [ "$1" = --remove ]; then
    rm -f "$dest"
    refresh --fresh
elif [ -s /run/secrets/ca ]; then
    cp /run/secrets/ca "$dest"
    refresh
fi
