#!/bin/bash

# compiles Jie Chen's RLCM package
cp -f configure_alexhack.sh RLCM
(cd RLCM; ./configure_alexhack.sh && make && make install)

# compile our wrapper around it
# *** to do
