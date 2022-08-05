#!/bin/bash

# compiles Jie Chen's RLCM package. We need the "apps", and we added
# wrappers to the SE kernel there which use IO via binary files...
cp -f configure_linux.sh RLCM
(cd RLCM; ./configure_linux.sh ${1:-g++} && make && make install && make apps)
