#!/bin/bash

# compiles Jie Chen's RLCM package. We need the "apps"
cp -f configure_alexhack.sh RLCM
(cd RLCM; ./configure_alexhack.sh && make && make install && make apps)

# compile our wrapper around it.
# *** need cpp wrapper which reads binary (train, test x), writes to binary
# (pred y, pred ytrg).  Since his .ex app files don't write out the answer!




