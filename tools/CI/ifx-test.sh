#!/usr/bin/env bash

source /opt/intel/oneapi/setvars.sh > /dev/null

folder="test-ifx"
rm -rf $folder && mkdir $folder && pushd $folder
FC=ifx cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON 
make
ctest --output-on-failure
