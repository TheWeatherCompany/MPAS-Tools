#!/bin/sh

module swap PrgEnv-cray PrgEnv-intel
module unload cray-netcdf
module load cray-netcdf-hdf5parallel
export NETCDF=$NETCDF_DIR
unlink mesh_to_mesh_interp
ln -sf Makefile.cray Makefile
make clean
make
if [ -f "a.out" ]; then
    ln -sf a.out ./mesh_to_mesh_interp
fi
