#include <iostream>
#include <typeinfo>
#include <typeindex>
#include <math.h>
#include "RemapperCell.h"
#include "array_utils.hpp"
#include "interp_utils.h"

#define HSrcPts2d(i,j) HSrcPts[(i)*maxHSrcPts+(j)]
#define HSrcWghts2d(i,j) HSrcWghts[(i)*maxHSrcPts+(j)]
#define nVSrcPts2d(i,j) nVSrcPts[(i)*nVDstPts+(j)]
#define VSrcPts3d(i,j,k) VSrcPts[((i)*nVDstPts+(j))*maxVSrcPts+(k)]
#define VSrcWghts3d(i,j,k) VSrcWghts[((i)*nVDstPts+(j))*maxVSrcPts+(k)]

RemapperCell::RemapperCell()
{
	nHDstPts = 0;
	maxHSrcPts = 0;
	nHSrcPts = NULL;
	HSrcPts = NULL;
	HSrcWghts = NULL;
	HSrcNearest = NULL;
	HSrcNearestLand = NULL;
	HSrcNearestWater = NULL;
	HSrcNearestSameLandmask = NULL;

	nVDstPts = 0;
	nVSrcLevels = 0;
	maxVSrcPts = 0;
	nVSrcPts = NULL;
	VSrcPts = NULL;
	VSrcWghts = NULL;
}

RemapperCell::RemapperCell(int nCellsDst, int nVertLevelsSrc, int nVertLevelsDst, int vertexDegree, int nSoilLevs)
{
    nHDstPts = nCellsDst;
    maxHSrcPts = vertexDegree;
    nHSrcPts = new int[nHDstPts];   // set to all 1 for now...
    HSrcPts = new int[(size_t)nHDstPts * (size_t)maxHSrcPts];
    HSrcWghts = new float[(size_t)nHDstPts * (size_t)maxHSrcPts];
    HSrcNearest = new int[(size_t)nHDstPts];
    HSrcNearestLand = new int[(size_t)nHDstPts];
    HSrcNearestWater = new int[(size_t)nHDstPts];
    HSrcNearestSameLandmask = new int[(size_t)nHDstPts];

    nVDstPts = nVertLevelsDst;
    nVSrcLevels = nVertLevelsSrc;
    maxVSrcPts = 2;
    nSoilLevels = nSoilLevs;
    
    if (nVSrcLevels > 0 && nVDstPts > 0) {
        nVSrcPts = new int[(size_t)nHDstPts * (size_t)nVDstPts];
        VSrcPts = new int[(size_t)nHDstPts * (size_t)nVDstPts * (size_t)maxVSrcPts];
        VSrcWghts = new float[(size_t)nHDstPts * (size_t)nVDstPts * (size_t)maxVSrcPts];
    }
}

RemapperCell::~RemapperCell()
{
	if (nHSrcPts != NULL) delete[] nHSrcPts;
	if (HSrcPts != NULL) delete[] HSrcPts;
	if (HSrcWghts != NULL) delete[] HSrcWghts;
	if (HSrcNearest != NULL) delete[] HSrcNearest;
	if (HSrcNearestLand != NULL) delete[] HSrcNearestLand;
	if (HSrcNearestWater != NULL) delete[] HSrcNearestWater;
	if (HSrcNearestSameLandmask != NULL) delete[] HSrcNearestSameLandmask;

	if (nVSrcPts != NULL) delete[] nVSrcPts;
	if (VSrcPts != NULL) delete[] VSrcPts;
	if (VSrcWghts != NULL) delete[] VSrcWghts;
}

void RemapperCell::computeWeightsCell(int nCellSrc, int *nEdgesOnCellSrc, int **verticesOnCellSrc, 
                                      int **cellsOnVertexSrc, int **cellsOnCellSrc,
                                      float *xCellSrc, float *yCellSrc, float *zCellSrc,
                                      float *xVertexSrc, float *yVertexSrc, float *zVertexSrc, float **levelsSrc,
                                      int *landmaskSrc,
                                      float *xCellDst, float *yCellDst, float *zCellDst, int *landmaskDst, 
                                      float **levelsDst)
{
	float tempLevels[nVSrcLevels];
	float vertCoords[maxHSrcPts][3];
	float pointInterp[3];
    bool do3d = (nVSrcLevels > 0 && nVDstPts > 0);

	int jCell = 0;
#pragma omp parallel for firstprivate(jCell) private(pointInterp, vertCoords, tempLevels)
	for (int i=0; i<nHDstPts; i++) {
		nHSrcPts[i] = maxHSrcPts;
        if (i % 1000000 == 0) {
            fprintf(stderr,"processing i: %d\n",i);
        }
		jCell = nearest_vertex(xCellDst[i], yCellDst[i], zCellDst[i], jCell, maxHSrcPts,
                           nEdgesOnCellSrc, verticesOnCellSrc, cellsOnVertexSrc,
                           xCellSrc, yCellSrc, zCellSrc, xVertexSrc, yVertexSrc, zVertexSrc);
		HSrcPts2d(i,0) = cellsOnVertexSrc[jCell][0] - 1;
		HSrcPts2d(i,1) = cellsOnVertexSrc[jCell][1] - 1;
		HSrcPts2d(i,2) = cellsOnVertexSrc[jCell][2] - 1;
        
        pointInterp[0] = xCellDst[i];
        pointInterp[1] = yCellDst[i];
        pointInterp[2] = zCellDst[i];
        for (int iv=0; iv<maxHSrcPts; iv++) {
            vertCoords[iv][0] = xCellSrc[HSrcPts2d(i,iv)];
            vertCoords[iv][1] = yCellSrc[HSrcPts2d(i,iv)];
            vertCoords[iv][2] = zCellSrc[HSrcPts2d(i,iv)];
        }
        mpas_wachspress_coordinates(maxHSrcPts, vertCoords, pointInterp, HSrcWghts + (i*maxHSrcPts));
        float maxWght = 0.0;
        float maxWghtLand = 0.0;
        float maxWghtWater = 0.0;
        float maxWghtSameLandmask = 0.0;
        int maxWghtPt = -1;
        int maxWghtPtLand = -1;
        int maxWghtPtWater = -1;
        int maxWghtPtSameLandmask = -1;
        for (int j = 0; j < nHSrcPts[i]; j++)
        {
            if (HSrcWghts2d(i,j) > maxWght)
            {
                maxWghtPt = HSrcPts2d(i,j);
                maxWght = *(HSrcWghts + (i*maxHSrcPts) + j);
            }
            if (landmaskSrc[HSrcPts2d(i,j)] == 1)
            {
                if (HSrcWghts2d(i,j) > maxWghtLand)
                {
                    maxWghtPtLand = HSrcPts2d(i,j);
                    maxWghtLand = *(HSrcWghts + (i*maxHSrcPts) + j);
                }
            }
            if (landmaskSrc[HSrcPts2d(i,j)] == 0)
            {
                if (HSrcWghts2d(i,j) > maxWghtWater)
                {
                    maxWghtPtWater = HSrcPts2d(i,j);
                    maxWghtWater = *(HSrcWghts + (i*maxHSrcPts) + j);
                }
            }
            if (landmaskSrc[HSrcPts2d(i,j)] == landmaskDst[i])
            {
                if (HSrcWghts2d(i,j) > maxWghtSameLandmask)
                {
                    maxWghtPtSameLandmask = HSrcPts2d(i,j);
                    maxWghtSameLandmask = *(HSrcWghts + (i*maxHSrcPts) + j);
                }
            }
        }

        HSrcNearest[i] = maxWghtPt;

        if (maxWghtPtLand != -1)
            HSrcNearestLand[i] = maxWghtPtLand;
        else
            HSrcNearestLand[i] = maxWghtPt;

        if (maxWghtPtWater != -1)
            HSrcNearestWater[i] = maxWghtPtWater;
        else
            HSrcNearestWater[i] = maxWghtPt;

        if (maxWghtPtSameLandmask != -1)
            HSrcNearestSameLandmask[i] = maxWghtPtSameLandmask;
        else
            HSrcNearestSameLandmask[i] = maxWghtPt;

		if (do3d) {
			// Horizontally interpolate column of levelsSrc values
			for (int k=0; k<nVSrcLevels; k++) {
				tempLevels[k] = 0;
			}
			for (int j=0; j<nHSrcPts[i]; j++) {
				for (int k=0; k<nVSrcLevels; k++) {
					tempLevels[k] += (HSrcWghts2d(i,j) * levelsSrc[HSrcPts2d(i,j)][k]);
				}
			}

			// For each vertical destination point, determine weights from tempLevels
			for (int k=0; k<nVDstPts; k++) {
				get_weights_1d(nVSrcLevels, tempLevels, levelsDst[i][k],
                                               nVSrcPts + (i*nVDstPts + k),
                                               VSrcPts + ((i*nVDstPts+k)*maxVSrcPts),
                                               VSrcWghts + ((i*nVDstPts+k)*maxVSrcPts));
			}
		}
	}
}


void RemapperCell::remap(const std::type_info& t, int ndims, interp_type interp, void *dst, void *src)
{
	if (std::type_index(t) == typeid(float)) {
		if (ndims == 1) {
			float *dstf = (float *)dst;
			float *srcf = (float *)src;
            if (interp == nearest)
            {
                remap1DNearest(dstf, srcf, HSrcNearest);
            } 
            else if (interp == nearest_land)
            {
                remap1DNearest(dstf, srcf, HSrcNearestLand);
            } 
            else if (interp == nearest_water)
            {
                remap1DNearest(dstf, srcf, HSrcNearestWater);
            } 
            else if (interp == nearest_samelandmask)
            {
                remap1DNearest(dstf, srcf, HSrcNearestSameLandmask);
            }
            else if (interp == barycentric) 
            {
                remap1D(dstf, srcf);
            } 
            else 
            {
                std::cout << "Unsupported interpolation type for 1D field: " << interp << std::endl;
            }
		}
		else if (ndims == 2) {
			float **dstf = (float **)dst;
			float **srcf = (float **)src;
            if (interp == nearest)
            {
                remap2DNearest(dstf, srcf, HSrcNearest);
            } 
            else if (interp == nearest_land)
            {
                remap2DNearest(dstf, srcf, HSrcNearestLand);
            } 
            else if (interp == nearest_water)
            {
                remap2DNearest(dstf, srcf, HSrcNearestWater);
            } 
            else if (interp == nearest_samelandmask)
            {
                remap2DNearest(dstf, srcf, HSrcNearestSameLandmask);
            }
            else if (interp == barycentric) 
            {
                remap2D(dstf, srcf);
            }
            else 
            {
                std::cout << "Unsupported interpolation type for 2D field: " << interp << std::endl;
            }
		}
		else if (ndims == 3) {
			float ***dstf = (float ***)dst;
			float ***srcf = (float ***)src;
            if (interp == barycentric)
            {
                remap3D(dstf, srcf);
            } 
            else if (interp == nearest_land)
            {
                remap3DNearest(dstf, srcf, HSrcNearestLand);
            }
            else if (interp == nearest_water)
            {
                remap3DNearest(dstf, srcf, HSrcNearestWater);
            }
            else if (interp == nearest_samelandmask)
            {
                remap3DNearest(dstf, srcf, HSrcNearestSameLandmask);
            }
            else if (interp == nearest_samelandmask_soil)
            {
                remap3DNearestSoil(dstf, srcf, HSrcNearestSameLandmask);
            }
            else 
            {
                std::cout << "Unsupported interpolation type for 3D field: " << interp << std::endl;
            }
		}
	}
	else {
		throw "RemapperCell can only handle 2-d or 3-d float fields";
	}
}

void RemapperCell::remap1D(float *dst, float *src)
{
	std::cerr << "Remapping 1d field\n";

#pragma omp parallel for
	for (int i=0; i<nHDstPts; i++) {
		dst[i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[i] += (HSrcWghts2d(i,j) * src[HSrcPts2d(i,j)]);
		}
	}
}

void RemapperCell::remap1DNearest(float *dst, float *src, int *nearestMap)
{
	std::cerr << "Remapping 1d field using nearest\n";

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for
    for (int i=0; i<nHDstPts; i++) {
        dst[i] = src[nearestMap[i]];
    }
}

void RemapperCell::remap2D(float **dst, float **src)
{
	std::cerr << "Remapping 2d field\n";

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for
	for (int i=0; i<nHDstPts; i++) {
		dst[0][i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[0][i] += (HSrcWghts2d(i,j) * src[0][HSrcPts2d(i,j)]);
		}
	}
}

void RemapperCell::remap2DNearest(float **dst, float **src, int *nearestMap)
{
	std::cerr << "Remapping 2d field using nearest\n";

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for
    for (int i=0; i<nHDstPts; i++) {
        dst[0][i] = src[0][nearestMap[i]];
    }
}


void RemapperCell::remap3D(float ***dst, float ***src)
{
	std::cerr << "Remapping 3d field\n";

	float tempLevels[nVSrcLevels];

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for private(tempLevels) schedule(dynamic,1000)
	for (size_t i=0; i<nHDstPts; i++) {
		// Horizontally interpolate column of levelsSrc values
		for (int k=0; k<nVSrcLevels; k++) {
			tempLevels[k] = 0;
		}
		for (int j=0; j<nHSrcPts[i]; j++) {
			for (int k=0; k<nVSrcLevels; k++) {
				tempLevels[k] += (HSrcWghts2d(i,j) * src[0][HSrcPts2d(i,j)][k]);
			}
		}

		// For each vertical destination point, interpolate
		for (size_t k=0; k<nVDstPts; k++) {
			dst[0][i][k] = 0;
			for (size_t j=0; j<nVSrcPts2d(i,k); j++) {
				dst[0][i][k] += VSrcWghts3d(i,k,j) * tempLevels[VSrcPts3d(i,k,j)];
			}
		}
	}
}

void RemapperCell::remap3DNearestSoil(float ***dst, float ***src, int *nearestMap)
{
    std::cerr << "Remapping 3d soil field for " << nSoilLevels << " levels. \n";

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for schedule(dynamic,1000)
	for (size_t i=0; i<nHDstPts; i++) {
		// Horizontally interpolate column of levelsSrc values
		for (int k=0; k<nSoilLevels; k++) {
			dst[0][i][k] = src[0][nearestMap[i]][k];
		}
	}
}

void RemapperCell::remap3DNearest(float ***dst, float ***src, int *nearestMap)
{
	std::cerr << "Remapping 3d field\n";

	float tempLevels[nVSrcLevels];

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for private(tempLevels) schedule(dynamic,1000)
	for (size_t i=0; i<nHDstPts; i++) {
		// Horizontally interpolate column of levelsSrc values
		for (int k=0; k<nVSrcLevels; k++) {
			tempLevels[k] = src[0][nearestMap[i]][k];
		}

		// For each vertical destination point, interpolate
		for (size_t k=0; k<nVDstPts; k++) {
			dst[0][i][k] = 0;
			for (size_t j=0; j<nVSrcPts2d(i,k); j++) {
				dst[0][i][k] += VSrcWghts3d(i,k,j) * tempLevels[VSrcPts3d(i,k,j)];
			}
		}
	}
}
