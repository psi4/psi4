set -eu -o pipefail

# Temporarily change directory to $HOME to install software
pushd .
cd "$HOME"

# Install Miniconda
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    # Make OSX md5 mimic md5sum from linux, alias does not work
    md5sum () {
        command md5 -r "$@"
    }
    MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh
else
    MINICONDA=Miniconda3-latest-Linux-x86_64.sh
fi
MINICONDA_HOME=$HOME/miniconda

echo "-- Installing latest Miniconda"
if [ -d "$MINICONDA_HOME/bin" ]; then
    echo "-- Miniconda latest version FOUND in cache"
    
    # Config settings are not saved in the cache
    export PIP_ARGS="-U"
    export PATH=$MINICONDA_HOME/bin:$PATH

    conda config --set always_yes yes --set changeps1 no
else
    MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
    echo "-- Miniconda latest version NOT FOUND in cache"
    wget -q https://repo.continuum.io/miniconda/$MINICONDA
    if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
        echo "Miniconda MD5 mismatch"
        exit 1
    fi
    # Travis creates the cached directories for us.
    # This is problematic when wanting to install Anaconda for the first time...
    rm -rf "$MINICONDA_HOME"
    bash $MINICONDA -b -p "$MINICONDA_HOME"

    # Configure miniconda
    export PIP_ARGS="-U"
    export PATH=$MINICONDA_HOME/bin:$PATH

    conda config --set always_yes yes --set changeps1 no
    conda update --q conda

    rm -f $MINICONDA
fi
echo "-- Done with latest Miniconda"

# Restore original directory
popd
