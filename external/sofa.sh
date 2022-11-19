#!/bin/bash

# Store current directory
DIR=$(pwd)

# Exit if SOFA already in directory
if test -f "libsofa.a"
then
    echo "SOFA already in directory"
    exit 0
fi

# Download and extract SOFA
curl -O http://www.iausofa.org/2021_0512_F/sofa_f-20210512.tar.Z
tar -zxvf sofa_f-20210512.tar.Z

# Compile SOFA
cd $DIR/sofa/20210512/f77/src
make

# Copy to working directory
cp libsofa.a $DIR

# Clear-up
rm -rf $DIR/sofa/
rm $DIR/sofa_f-20210512.tar.Z