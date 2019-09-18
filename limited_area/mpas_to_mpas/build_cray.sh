#!/bin/sh

module swap PrgEnv-cray PrgEnv-intel
module unload cray-netcdf
module load cray-netcdf-hdf5parallel
export NETCDF=$NETCDF_DIR

# mesh_to_mesh_interp
echo "Building mesh_to_mesh_interp"
unlink mesh_to_mesh_interp
ln -sf Makefile.cray Makefile
make clean
make
if [ -f "a.out" ]; then
    cp a.out ./mesh_to_mesh_interp
fi

# merge_fields
echo "Building merge_fields"
ftn merge_interpolated_variables.f90 -L${NETCDF}/lib -o merge_fields
