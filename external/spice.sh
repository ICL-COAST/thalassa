#!/bin/bash

# Store current directory
DIR=$(pwd)

# Check if SPICE already present
echo "Checking for SPICE..."

# Exit if SPICE already in directory
if test -f "spicelib.a"
then
    echo "- Library already in directory: spicelib.a"
    exit 0
fi

# Download and extract SPICE
echo "Downloading SPICE..."
curl -O https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/PC_Linux_gfortran_64bit/packages/toolkit.tar.Z > /dev/null 2>&1
tar -zxvf toolkit.tar.Z > /dev/null 2>&1

# Compile SPICE
echo "Compiling SPICE..."
cd $DIR/toolkit
./makeall.csh > /dev/null 2>&1

# Copy to working directory
cp $DIR/toolkit/lib/spicelib.a $DIR

# Clear-up
echo "Cleaning up..."
rm -rf $DIR/toolkit/
rm $DIR/toolkit.tar.Z