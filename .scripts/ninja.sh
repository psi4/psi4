#!/usr/bin/env bash

set -eu -o pipefail

echo "-- Installing Ninja"
Ninja_VERSION="1.9.0.g99df1.kitware.dyndep-1.jobserver-1"
Ninja_URL="https://github.com/Kitware/ninja/releases/download/v${Ninja_VERSION}/ninja-${Ninja_VERSION}_x86_64-linux-gnu.tar.gz"
cd "$HOME"/Deps
mkdir -p ninja
curl -Ls $Ninja_URL | tar -xz -C ninja --strip-components=1
cd "$TRAVIS_BUILD_DIR"
echo "-- Done with Ninja"
