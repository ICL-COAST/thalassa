#!/bin/bash

# Store current directory
DIR=$(pwd)

# Exit if SPICE already in directory
if test -f "spicelib.a"
then
    echo "SPICE already in directory"
    exit 0
fi

# Download and extract SPICE
curl -O https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/PC_Linux_gfortran_64bit/packages/toolkit.tar.Z
tar -zxvf toolkit.tar.Z

# Compile SPICE
cd $DIR/toolkit
./makeall.csh

# Copy to working directory
cp $DIR/toolkit/lib/spicelib.a $DIR

# Clear-up
rm -rf $DIR/toolkit/
rm $DIR/toolkit.tar.Z