#!/bin/bash

# compiles Jie Chen's RLCM package. We need the "apps", and we added
# wrappers to the SE kernel there which use IO via binary files...
cp -f configure_alexhack.sh RLCM
(cd RLCM; ./configure_alexhack.sh && make && make install && make apps)
