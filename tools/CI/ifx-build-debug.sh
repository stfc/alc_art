#!/usr/bin/env bash

source /opt/intel/oneapi/setvars.sh > /dev/null

folder="build-ifx-debug"
rm -rf $folder && mkdir $folder && pushd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
