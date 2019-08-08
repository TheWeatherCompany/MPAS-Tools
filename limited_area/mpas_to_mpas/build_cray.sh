#!/bin/sh

module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf
export NETCDF=$NETCDF_DIR
unlink mesh_to_mesh_interp
make clean
make
if [ -f "a.out" ]; then
    ln -sf a.out ./mesh_to_mesh_interp
fi
