#!/bin/bash

if [ "$#" -lt 3 ]; then
   echo "USAGE: ./run_pbsim.sh ref.fa sample.fq depth"
   exit
fi

getFullDir_file(){
    echo "$(cd "$(dirname "$1")"; pwd)"
}

SCRIPT_DIR=`getFullDir_file "${BASH_SOURCE[0]}"`

${SCRIPT_DIR}/pbsim/src/pbsim \
--data-type CLR \
--depth ${3} \
--length-min 1 \
--length-max 100000 \
--seed 0 \
--sample-fastq ${2} \
${1} \
