#!/bin/bash
 
# Declare kernel paths
kernelarray=(
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_part-1.bsp"
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de431_part-2.bsp"
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc"
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
)

# Declare empty string
downloadstr=""

echo "Checking for downloaded kernels..."

# Iterate through kernels
for kernel in ${kernelarray[*]}; do
    # Strip the filename
    kernelbasename="${kernel##*/}"

    # Check if kernel already downloaded
    if test -f $kernelbasename
    then
        echo "- Kernel already in directory: $kernelbasename"
    else
        echo "- Kernel to be downloaded:     $kernelbasename"
        downloadstr="$downloadstr -O $kernel"
    fi
done

# Download kernels
if [ -z "${downloadstr}" ]
then
    exit 0
else
    echo "Downloading remaining kernels..."
    curl -Z $downloadstr
fi