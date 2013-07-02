#!/bin/sh

# Warning: mkdir -p isn't portable

if [ ! -d $* ]; then
  mkdir -p $*
fi
