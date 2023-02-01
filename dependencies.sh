#!/bin/bash

pwd=$(pwd)

# Download kernels
cd $pwd/data/kernels/
./kernels.sh
echo ""