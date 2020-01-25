#!/usr/bin/env python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description='Assign MPAS cell points to convex partitions specified by centroids.')

parser.add_argument('-g', action='store', dest='grid_file', default='grid.nc',
                    help='NetCDF file containing {x,y,z}Cell. (default: grid.nc)')
parser.add_argument('centroids', action='store', nargs='+',
                    help='Text file(s) containing cartesian coordinates for partitioning centroids')

arguments = parser.parse_args()

grid_file = arguments.grid_file
centroid_files = arguments.centroids

from netCDF4 import Dataset as NetCDFFile
from scipy.spatial import cKDTree as KDTree
import numpy as np

with NetCDFFile(grid_file) as nc:
    xCell = nc.variables['xCell'][:]
    yCell = nc.variables['yCell'][:]
    zCell = nc.variables['zCell'][:]
cells = np.vstack([xCell,yCell,zCell]).T
print("Number of cells: {:d}".format(len(xCell)))

for i, centroid_file in enumerate(centroid_files):
    centroids = np.loadtxt(centroid_file)
    print("[{:d}/{:d}] Number of centroids in {}: {:d}".format(i+1, len(centroid_files), centroid_file, len(centroids)))
    _, nn = KDTree(centroids).query(cells)
    nn = nn + 1 # Switch to Fortran indexing
    _, counts = np.unique(nn, return_counts=True)
    print("        Smallest partition: {:d}".format(np.min(counts)))
    print("        Largest partition: {:d} (+{:d})".format(np.max(counts), np.max(counts)-np.min(counts)))
    dest = centroid_file + '.cvt'
    print("        Saving to " + dest)
    np.savetxt(dest, nn, fmt='%d')

print("Done!")
