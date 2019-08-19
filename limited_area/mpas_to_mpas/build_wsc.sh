#!/bin/sh

unlink mesh_to_mesh_interp
ln -sf Makefile.wsc Makefile
make clean
make
if [ -f "a.out" ]; then
    ln -sf a.out ./mesh_to_mesh_interp
fi
