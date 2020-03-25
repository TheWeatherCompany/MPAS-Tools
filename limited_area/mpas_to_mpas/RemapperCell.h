#ifndef _REMAPPERCELL_H
#define _REMAPPERCELL_H

#include "RemapperBase.h"

class RemapperCell : virtual public RemapperBase {
public:
	RemapperCell();
    RemapperCell(int nCellsDst, int nVertLevelsSrc, int nVertLevelsDst, int vertexDegree);
	~RemapperCell();
	void remap(const std::type_info& t, int ndims, void *dst, void *src);
	void computeWeightsCell(int *nEdgesOnCellSrc, int **verticesOnCellSrc, int **cellsOnVertexSrc,
                            float *xCellSrc, float *yCellSrc, float *zCellSrc,
                            float *xVertexSrc, float *yVertexSrc, float *zVertexSrc,
                            float **levelsSrc,
                            float *xCellDst, float *yCellDst, float *zCellDst,
                            float **levelsDst);

private:
	//
	// Horizontal remapping fields
	//
	size_t nHDstPts;      // Number of horizontal destination points
	size_t maxHSrcPts;    // Maximum number of horizontal source points needed by any destination point
	int *nHSrcPts;     // Number of horizontal source points needed by any destination point
	int *HSrcPts;      // Source points needed for horizontal interpolation to each destination point
	float *HSrcWghts;  // Source weights needed for horizontal interpolation to each destination point
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
	void remap2D(float **dst, float **src);
	void remap3D(float ***dst, float ***src);
};
#endif
