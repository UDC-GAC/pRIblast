#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: preprocess.sh /path/to/dataset"
    exit 1
fi

DATASET_PATH=$(dirname "$1")
DATASET_NAME=$(basename "$1")

cat $1 |
    while read line
    do
        if [ "${line:0:1}" == ">" ]
        then
            echo -e "\n"$line
        else
            echo $line | tr -d '\n'
        fi
    done |
        tail -n+2 > db-$DATASET_NAME && grep -A 1 lncRNA db-$DATASET_NAME > ris-$DATASET_NAME
