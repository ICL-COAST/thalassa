 # Setup of HDF5
 HDF5 is a data format that is used by THALASSA to save data in batch mode. To properly run THALASSA with HDF5, one first needs to install the HDF5 library.

 ## Installing the HDF5 library
 For MacOS, download the [source](https://www.hdfgroup.org/downloads/hdf5/source-code/#conf "HDF5 Downloads") with Unix line endings. Decompress the archive in a directory of your choice, for instance /Users/davide/hdf5-1.XX.Y. If you want to install the HDF5 library in the directory /Users/davide/hdf5/, navigate to the directory and issue the command:
      
    ./configure --prefix=/Users/davide/hdf5/ --enable-fortran

Once this is done, run:

    make
    make check  # have a cup of coffee because this will take a while
    make install
    make check-install

