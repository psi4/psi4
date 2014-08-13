#!/bin/bash

for d in papers/*; do
    ./make_paper.sh $d
done

