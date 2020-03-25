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
#deinfe VSrcPts3d(i,j,k) VSrcPts[(i)*nVDstPts*maxVSrcPts+(j)*maxVSrcPts+(k)]
#deinfe VSrcWghts3d(i,j,k) VSrcWghts[(i)*nVDstPts*maxVSrcPts+(j)*maxVSrcPts+(k)]

RemapperCell::RemapperCell()
{
	nHDstPts = 0;
	maxHSrcPts = 0;
	nHSrcPts = NULL;
	HSrcPts = NULL;
	HSrcPts2d = NULL;
	HSrcWghts = NULL;
	HSrcWghts2d = NULL;

	nVDstPts = 0;
	nVSrcLevels = 0;
	maxVSrcPts = 0;
	nVSrcPts = NULL;
	nVSrcPts2d = NULL;
	VSrcPts = NULL;
	VSrcPts3d = NULL;
	VSrcWghts = NULL;
	VSrcWghts3d = NULL;
}

RemapperCell::RemapperCell(int nCellsDst, int nVertLevelsSrc, int nVertLevelsDst, int vertexDegree)
{
    nHDstPts = nCellsDst;
    maxHSrcPts = vertexDegree;
    nHSrcPts = new int[nHDstPts];   // set to all 1 for now...
    HSrcPts = new int[(size_t)nHDstPts * (size_t)maxHSrcPts];
    HSrcWghts = new float[(size_t)nHDstPts * (size_t)maxHSrcPts];
    
    nVDstPts = nVertLevelsDst;
    nVSrcLevels = nVertLevelsSrc;
    maxVSrcPts = 2;
    
    if (nVSrcLevels > 0 && nVDstPts > 0) {
        nVSrcPts = new int[(size_t)nHDstPts * (size_t)nVDstPts];
        nVSrcPts2d = allocate_2d<int>(nHDstPts, nVDstPts, nVSrcPts);
        VSrcPts = new int[(size_t)nHDstPts * (size_t)nVDstPts * (size_t)maxVSrcPts];
        VSrcWghts = new float[(size_t)nHDstPts * (size_t)nVDstPts * (size_t)maxVSrcPts];
    }
}

RemapperCell::~RemapperCell()
{
	if (nHSrcPts != NULL) delete[] nHSrcPts;
	if (HSrcPts != NULL) delete[] HSrcPts;
	if (HSrcWghts != NULL) delete[] HSrcWghts;

	if (nVSrcPts != NULL) delete[] nVSrcPts;
	if (VSrcPts != NULL) delete[] VSrcPts;
	if (VSrcWghts != NULL) delete[] VSrcWghts;
}

void RemapperCell::computeWeightsCell(int *nEdgesOnCellSrc, int **verticesOnCellSrc, int **cellsOnVertexSrc,
                                      float *xCellSrc, float *yCellSrc, float *zCellSrc,
                                      float *xVertexSrc, float *yVertexSrc, float *zVertexSrc, float **levelsSrc,
                                      float *xCellDst, float *yCellDst, float *zCellDst, float **levelsDst)
{
	int j;
	int nCells;
	float tempLevels[nVSrcLevels];
	float vertCoords[maxHSrcPts][3];
	float pointInterp[3];
    bool do3d = (nVSrcLevels > 0 && nVDstPts > 0);

	int jCell = 0;
#pragma omp parallel for firstprivate(jCell) private(pointInterp, vertCoords, tempLevels)
	for (int i=0; i<nHDstPts; i++) {
		nHSrcPts[i] = maxHSrcPts;
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
            vertCoords[iv][0] = xCellSrc[HSrcPts2d[i][iv]];
            vertCoords[iv][1] = yCellSrc[HSrcPts2d[i][iv]];
            vertCoords[iv][2] = zCellSrc[HSrcPts2d[i][iv]];
        }

		mpas_wachspress_coordinates(maxHSrcPts, vertCoords, pointInterp, HSrcWghts2d[i]);

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
				get_weights_1d(nVSrcLevels, tempLevels, levelsDst[i][k], &nVSrcPts2d(i,k), &VSrcPts3d(i,k,0), &VSrcWghts3d(i,k,0));
			}
		}
	}
}


void RemapperCell::remap(const std::type_info& t, int ndims, void *dst, void *src)
{
	if (std::type_index(t) == typeid(float)) {
		if (ndims == 1) {
			float *dstf = (float *)dst;
			float *srcf = (float *)src;
			remap1D(dstf, srcf);
		}
		else if (ndims == 2) {
			float **dstf = (float **)dst;
			float **srcf = (float **)src;
			remap2D(dstf, srcf);
		}
		else if (ndims == 3) {
			float ***dstf = (float ***)dst;
			float ***srcf = (float ***)src;
			remap3D(dstf, srcf);
		}
	}
	else {
		throw "RemapperCell can only handle 2-d or 3-d float fields";
	}
}

void RemapperCell::remap1D(float *dst, float *src)
{
	std::cerr << "Remapping 1d field\n";

	for (int i=0; i<nHDstPts; i++) {
		dst[i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[i] += (HSrcWghts2d(i,j) * src[HSrcPts2d(i,j)]);
		}
	}
}

void RemapperCell::remap2D(float **dst, float **src)
{
	std::cerr << "Remapping 2d field\n";

	// TODO: Right now, the time dimension is the first dimension
	for (int i=0; i<nHDstPts; i++) {
		dst[0][i] = 0;
		for (int j=0; j<nHSrcPts[i]; j++) {
			dst[0][i] += (HSrcWghts2d(i,j) * src[0][HSrcPts2d(i,j)]);
		}
	}
}

void RemapperCell::remap3D(float ***dst, float ***src)
{
	std::cerr << "Remapping 3d field\n";

	float tempLevels[nVSrcLevels];

	// TODO: Right now, the time dimension is the first dimension
#pragma omp parallel for private(tempLevels) schedule(dynamic,1000)
	for (int i=0; i<nHDstPts; i++) {
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
		for (int k=0; k<nVDstPts; k++) {
			dst[0][i][k] = 0;
			for (int j=0; j<nVSrcPts2d(i,k); j++) {
				dst[0][i][k] += VSrcWghts3d(i,k,j) * tempLevels[VSrcPts3d(i,k,j)];
			}
		}
	}
}
