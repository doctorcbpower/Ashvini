#!/bin/bash
# Make the two tarballs to copy to a cluster without internet: Ashvini at HEAD and foraois at the commit in hpc/foraois_commit.txt.
# Usage: bash hpc/make_package.sh <path-to-foraois>
set -euo pipefail
FORAOIS_DIR="${1:?usage: bash hpc/make_package.sh <path-to-foraois>}"
COMMIT=$(cut -d' ' -f1 hpc/foraois_commit.txt)
git archive --prefix=Ashvini/ -o "ashvini_$(git rev-parse --short HEAD).tar.gz" HEAD
git -C "$FORAOIS_DIR" archive --prefix=foraois/ -o "$PWD/foraois_${COMMIT:0:7}.tar.gz" "$COMMIT"
sha256sum ashvini_*.tar.gz foraois_*.tar.gz > SHA256SUMS
echo "copy the two tarballs and SHA256SUMS; check with: sha256sum -c SHA256SUMS"
