#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: preprocess.sh dataset-name"
    exit 1
fi

DATASET=$1

cat $DATASET |
    while read line
    do
        if [ "${line:0:1}" == ">" ]
        then
            echo -e "\n"$line
        else
            echo $line | tr -d '\n'
        fi
    done |
        tail -n+2 > db-$DATASET && grep -A 1 lncRNA db-$DATASET > ris-$DATASET
