#!/bin/sh

# mesh_to_mesh_interp
echo "Building mesh_to_mesh_interp"
unlink mesh_to_mesh_interp
ln -sf Makefile.wsc Makefile
make clean
make
if [ -f "a.out" ]; then
    ln -sf a.out ./mesh_to_mesh_interp
fi

# merge_fields
echo "Building merge_fields"
mpif90 merge_interpolated_variables.f90 -L${NETCDF}/lib -o merge_fields
