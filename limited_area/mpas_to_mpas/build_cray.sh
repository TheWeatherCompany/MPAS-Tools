#!/bin/sh

module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf
export NETCDF=$NETCDF_DIR
unlink mesh_to_mesh_interp
ln -sf Makefile.cray Makefile
make clean
make
if [ -f "a.out" ]; then
    ln -sf a.out ./mesh_to_mesh_interp
fi
