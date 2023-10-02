#!/bin/bash
###############################################
### load your build environment 	    ###
###############################################
#----------------------------------------------------------------------------------#
# Setting for compilers
echo "Module Load"

module purge

if [[ `uname -n` == "chama"* ]] || [[ `uname -n` == "skybridge"* ]] || [[ `uname -n` == "cee-build"* ]]; then

  module load cde/v3/cmake/3.23.1
  module load cde/v3/intel-oneapi-compilers/2021.1.2
  module load cde/v3/intel-oneapi-mkl/2021.4.0-intel-2021.1.2
  module load cde/v3/intel-oneapi-mpi/2021.4.0-intel-2021.1.2
  module load cde/v3/hdf5/1.10.6-intel-2021.1.2-intel-oneapi-mpi-2021.4.0
  module load cde/v3/netcdf-c/4.8.1-intel-2021.1.2-intel-oneapi-mpi-2021.4.0
  module load cde/v3/boost/1.79.0-intel-2021.1.2-intel-oneapi-mpi-2021.4.0
  module load cde/v3/metis/5.1.0-intel-2021.1.2
  module load cde/v3/parmetis/4.0.3-intel-2021.1.2-intel-oneapi-mpi-2021.4.0

  export NETPUB=/projects/netpub
else

  #export NETPUB=/usr/netpub
  echo "This script only works for chama"
  exit 1
fi

DEPDIR=$NETPUB/Delft3D/Dependencies
#DELFTDIR=/projects/netpub/Delft3D/SNL-Delft3D-FM-CEC
#PWD=`pwd`

## not sure what this does, but seems to be required for mpiexec on chama
## credit: https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/Intel-MPI-Unable-to-run-bstrap-proxy-error-setting-up-the/m-p/1204677
export I_MPI_HYDRA_IFACE="ib0"

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
export F77=$FC
export F90=$FC


#----------------------------------------------------------------------------------#
# Environment variables for required packages

## patchelf
export PATCHELF=$DEPDIR/patchelf

export NETCDFFORT=$DEPDIR/netcdf-fortran-install
export PETSCPATH=$DEPDIR/petsc-install
export SQL3DIR=$DEPDIR/sqlite-install
export PROJDIR=$DEPDIR/proj-install
export GDALDIR=$DEPDIR/gdal-install

export PATH=$PATCHELF/bin:$PATH


export PKG_CONFIG_PATH=$NETCDFFORT/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$NETCDFFORT/lib:$LD_LIBRARY_PATH

export PKG_CONFIG_PATH=$PETSCPATH/lib/pkgconfig:$PKG_CONFIG_PATH

export PKG_CONFIG_PATH=$SQL3DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$SQL3DIR/lib:$LD_LIBRARY_PATH

export PKG_CONFIG_PATH=$PROJDIR/lib64/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$PROJDIR/lib64:$LD_LIBRARY_PATH

export PKG_CONFIG_PATH=$GDALDIR/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$GDALDIR/lib:$LD_LIBRARY_PATH

#echo $LD_LIBRARY_PATH
