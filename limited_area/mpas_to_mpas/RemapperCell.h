#ifndef _REMAPPERCELL_H
#define _REMAPPERCELL_H

#include "RemapperBase.h"

class RemapperCell : virtual public RemapperBase {
public:
	RemapperCell();
    RemapperCell(int nCellsDst, int nVertLevelsSrc, int nVertLevelsDst, int vertexDegree);
	~RemapperCell();
	void remap(const std::type_info& t, int ndims, interp_type interp, void *dst, void *src);
	void computeWeightsCell(int nCellSrc, int *nEdgesOnCellSrc, int **verticesOnCellSrc, 
                            int **cellsOnVertexSrc,
                            int **cellsOnCellSrc,
                            float *xCellSrc, float *yCellSrc, float *zCellSrc,
                            float *xVertexSrc, float *yVertexSrc, float *zVertexSrc,
                            float **levelsSrc, int *landmaskSrc,
                            float *xCellDst, float *yCellDst, float *zCellDst,
                            int *landmaskDst, float **levelsDst );

private:
	//
	// Horizontal remapping fields
	//
	size_t nHDstPts;      // Number of horizontal destination points
	size_t maxHSrcPts;    // Maximum number of horizontal source points needed by any destination point
	int *nHSrcPts;     // Number of horizontal source points needed by any destination point
	int *HSrcPts;      // Source points needed for horizontal interpolation to each destination point
	float *HSrcWghts;  // Source weights needed for horizontal interpolation to each destination point
    int *HSrcNearest; // Nearest source point to each destination point, used for nearest interp
    int *HSrcNearestLand; // Nearest source land point from surrounding trio of point to each destination point, if no land point exists this will be nearest point
    int *HSrcNearestWater; // Nearest source water point from surrounding trio of point to each destination point, if no water point exists this will be nearest point
    int *HSrcNearestSameLandmask; // Nearest source point from surrounding trio of point to each destination point that has the same landmask as the destination point. If there are not points with same mask, nearest point
	unsigned char *HDstMask;     // Mask for horizontal destination points

	//
	// Vertical remapping fields
	//
	size_t nVDstPts;      // Number of vertical destination points
	size_t nVSrcLevels;   // Number of vertical source points
	size_t maxVSrcPts;    // Maximum number of vertical source points needed by any destination point
	int *nVSrcPts;     // Number of vertical source points needed by any destination point
	int *VSrcPts;      // Source points needed for vertical interpolation to each destination point
	float *VSrcWghts;  // Source weights needed for vertical interpolation to each destination point

	//
	// Internal for handling remapping of 2- and 3-d fields
	//
	void remap1D(float *dst, float *src);
	void remap1DNearest(float *dst, float *src, int *nearestMap);
	void remap2D(float **dst, float **src);
	void remap2DNearest(float **dst, float **src, int *nearestMap);
	void remap3D(float ***dst, float ***src);
};
#endif
