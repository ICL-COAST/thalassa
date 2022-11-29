#!/bin/bash

pwd=$(pwd)

# Download SOFA
cd $pwd/external
./sofa.sh
echo ""

# Download SPICE
cd $pwd/external
./spice.sh
echo ""

# Download kernels
cd $pwd/data/kernels/
./kernels.sh
echo ""