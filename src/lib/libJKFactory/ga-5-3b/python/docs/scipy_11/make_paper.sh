#!/bin/bash

DIR=$1
WD=`pwd`

if [[ ! -d $DIR ]]; then
  echo "Usage: make_paper.sh source_dir"
  exit -1
fi

if [[ ! -d output ]]; then
  mkdir output
fi

AUTHOR=`basename $DIR`
OUTDIR="output/$AUTHOR"

mkdir -p $OUTDIR
cp $DIR/* $OUTDIR
python publisher/build_paper.py $DIR $OUTDIR
if [ "$?" -ne "0" ]; then
    echo "Error building paper $DIR. Aborting."
    exit 1
fi
cd $OUTDIR
pdflatex paper.tex && pdflatex paper.tex | (python $WD/publisher/paper_cut.py)

