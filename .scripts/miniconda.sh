#!/usr/bin/env bash

set -eu -o pipefail

echo "-- Installing latest Miniconda"
if [ -d "$HOME/Deps/miniconda/bin" ]; then
    echo "-- Miniconda latest version FOUND in cache"
else
    cd "$HOME"/Downloads
    echo "-- Miniconda latest version NOT FOUND in cache"
    if [ "$TRAVIS_OS_NAME" = "linux" ]; then
        curl -Ls https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh
    elif [ "$TRAVIS_OS_NAME" = "osx" ]; then
        curl -Ls https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh > miniconda.sh
    fi
    # Travis creates the cached directories for us.
    # This is problematic when wanting to install Anaconda for the first time...
    rm -rf "$HOME"/Deps/miniconda
    bash miniconda.sh -b -p "$HOME"/Deps/miniconda > /dev/null 2>&1
    touch "$HOME"/Deps/miniconda/conda-meta/pinned
    PATH=$HOME/Deps/miniconda/bin${PATH:+:$PATH}
    hash -r
    conda config --set show_channel_urls True > /dev/null 2>&1
    conda config --set always_yes yes > /dev/null 2>&1
    conda config --set changeps1 no > /dev/null 2>&1
    conda update --all --yes > /dev/null 2>&1
    conda clean -tipy > /dev/null 2>&1
    conda info -a
    cd "$TRAVIS_BUILD_DIR"
    rm -f "$HOME"/Downloads/miniconda.sh
fi
echo "-- Done with latest Miniconda"
