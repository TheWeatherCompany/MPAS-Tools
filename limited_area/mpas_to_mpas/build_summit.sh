#!/bin/sh

# mesh_to_mesh_interp
echo "Building mesh_to_mesh_interp"
unlink mesh_to_mesh_interp
ln -sf Makefile.summit Makefile
make clean
make
if [ -f "a.out" ]; then
    cp a.out ./mesh_to_mesh_interp
fi

