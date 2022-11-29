#!/bin/bash

# Store current directory
DIR=$(pwd)

# Check if SOFA already present
echo "Checking for SOFA..."

# Exit if SOFA already in directory
if test -f "libsofa.a"
then
    echo "- Library already in directory: libsofa.a"
    exit 0
fi

# Download and extract SOFA
echo "Downloading SOFA..."
curl -O http://www.iausofa.org/2021_0512_F/sofa_f-20210512.tar.Z > /dev/null 2>&1
tar -zxvf sofa_f-20210512.tar.Z > /dev/null 2>&1

# Compile SOFA
echo "Compiling SOFA..."
cd $DIR/sofa/20210512/f77/src
make > /dev/null 2>&1

# Copy to working directory
cp libsofa.a $DIR

# Clear-up
echo "Cleaning up..."
rm -rf $DIR/sofa/
rm $DIR/sofa_f-20210512.tar.Z