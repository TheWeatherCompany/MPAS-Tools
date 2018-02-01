#include <iostream>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include "NCField.hpp"
#include "RemapperCell.h"
#include "RemapperEdge.h"
#include "mpas_utils.h"

void start_timer(int n);
void stop_timer(int n, int *secs, int *n_secs);


int main(int argc, char **argv)
{
	int ncid;
	int stat;

	const int NUM_SCALARS = 6;
	char qxNames[NUM_SCALARS][3] = {"qv", "qc", "qr", "qi", "qs", "qg"};

	NCField<float> *latCellDst;
	NCField<float> *lonCellDst;
	NCField<float> *latEdgeDst;
	NCField<float> *lonEdgeDst;
	NCField<float> *angleEdgeDst;
	NCField<int> *cellsOnEdgeDst;
	NCField<int> *bdyMaskCellDst;
	NCField<int> *bdyMaskEdgeDst;
	NCField<float> *zgridDst;
	NCField<float> *zmidDst;
	NCField<float> *zedgeDst;
	NCField<float> *uDst;
	NCField<float> *vDst;
	NCField<float> *thetaDst;
	NCField<float> *rhoDst;
	NCField<float> *wDst;
	NCField<float> *qxDst[NUM_SCALARS];
	NCField<char> *xtime;
	float ***uDstArr;
	float ***vDstArr;
	float **zgridDstArr;
	float **zmidDstArr;
	float ***rhoDstArr;
	float *angleEdgeDstArr;

	char qxLbcName[64];

	const char *globalMeshFile;
	const char *regionalMeshFile;
	const char *globalFieldFile;
	char regionalFieldFile[64];
	char **xtimeArr;
	char date[14];

	NCField<float> *latCellSrc;
	NCField<float> *lonCellSrc;
	NCField<float> *latEdgeSrc;
	NCField<float> *lonEdgeSrc;
	NCField<float> *latVertexSrc;
	NCField<float> *lonVertexSrc;
	NCField<int> *nEdgesOnCellSrc;
	NCField<int> *nEdgesOnEdgeSrc;
	NCField<int> *cellsOnCellSrc;
	NCField<int> *verticesOnCellSrc;
	NCField<int> *cellsOnVertexSrc;
	NCField<int> *cellsOnEdgeSrc;
	NCField<int> *edgesOnCellSrc;
	NCField<int> *edgesOnEdgeSrc;
	NCField<float> *weightsOnEdgeSrc;
	NCField<float> *angleEdgeSrc;
	NCField<float> *zgridSrc;
	NCField<float> *zmidSrc;
	NCField<float> *zedgeSrc;
	NCField<float> *uSrc;
	NCField<float> *vSrc;
	NCField<float> *thetaSrc;
	NCField<float> *rhoSrc;
	NCField<float> *wSrc;
	NCField<float> *qxSrc[NUM_SCALARS];
	float ***uSrcArr;
	float ***vSrcArr;
	float **zgridSrcArr;
	float **zmidSrcArr;
	float ***rhoSrcArr;
	int *nEdgesOnEdgeSrcArr;
	int **edgesOnEdgeSrcArr;
	float **weightsOnEdgeSrcArr;
	float *angleEdgeSrcArr;
	RemapperCell *cellLayerMap;
	RemapperCell *cellLevelMap;
	RemapperCell *cellToEdgeMap;
	RemapperEdge *edgeMap;
	int secs, nsecs;
	int itime;
	int argv_idx;
	int use_reconstruct_winds;


	if (argc < 4) {
		std::cerr << "\nUsage: " << argv[0] << " [--use-reconstruct-winds] <global_IC_file> <regional_IC_file> <global_fields_file>+\n\n";
		return 1;
	}

	if (strcmp(argv[1], "--use-reconstruct-winds") == 0) {
		std::cout << "Using reconstructed winds at cell centers...\n";
		argv_idx = 2;
		use_reconstruct_winds = 1;
	}
	else {
		//
		// Try to catch other uses options besides --use-reconstruct-winds
		//
		if (argv[1][0] == '-') {
			std::cerr << "Unrecognized option: " << argv[1] << std::endl;
			std::cerr << "\nSupported options are: --use-reconstruct-winds\n";
			return 1;
		}

		std::cout << "Using normal component of winds at edges...\n";
		argv_idx = 1;
		use_reconstruct_winds = 0;
	}
	globalMeshFile = argv[argv_idx++];
	regionalMeshFile = argv[argv_idx++];


	//
	// Read mesh description fields from regional IC file
	//
	start_timer(0);
	try {
		latCellDst = new NCField<float>(regionalMeshFile, "latCell");
	}
	catch (int e) {
		std::cerr << "Error reading latCell field from " << globalMeshFile << std::endl;
		return 1;
	}
	lonCellDst = new NCField<float>(regionalMeshFile, "lonCell");
	latEdgeDst = new NCField<float>(regionalMeshFile, "latEdge");
	lonEdgeDst = new NCField<float>(regionalMeshFile, "lonEdge");
	angleEdgeDst = new NCField<float>(regionalMeshFile, "angleEdge");
	cellsOnEdgeDst = new NCField<int>(regionalMeshFile, "cellsOnEdge");
	bdyMaskCellDst = new NCField<int>(regionalMeshFile, "bdyMaskCell");
	bdyMaskEdgeDst = new NCField<int>(regionalMeshFile, "bdyMaskEdge");
	zgridDst = new NCField<float>(regionalMeshFile, "zgrid");
	stop_timer(0, &secs, &nsecs);
	printf("Time to read mesh fields from %s : %i.%9.9i\n", regionalMeshFile, secs, nsecs);


	//
	// Read mesh description fields from global IC file
	//
	start_timer(0);
	try {
		latCellSrc = new NCField<float>(globalMeshFile, "latCell");
	}
	catch (int e) {
		std::cerr << "Error reading latCell field from " << globalMeshFile << std::endl;
		return 1;
	}
	lonCellSrc = new NCField<float>(globalMeshFile, "lonCell");
	latEdgeSrc = new NCField<float>(globalMeshFile, "latEdge");
	lonEdgeSrc = new NCField<float>(globalMeshFile, "lonEdge");
	latVertexSrc = new NCField<float>(globalMeshFile, "latVertex");
	lonVertexSrc = new NCField<float>(globalMeshFile, "lonVertex");
	nEdgesOnCellSrc = new NCField<int>(globalMeshFile, "nEdgesOnCell");
	nEdgesOnEdgeSrc = new NCField<int>(globalMeshFile, "nEdgesOnEdge");
	weightsOnEdgeSrc = new NCField<float>(globalMeshFile, "weightsOnEdge");
	edgesOnEdgeSrc = new NCField<int>(globalMeshFile, "edgesOnEdge");
	angleEdgeSrc = new NCField<float>(globalMeshFile, "angleEdge");
	cellsOnCellSrc = new NCField<int>(globalMeshFile, "cellsOnCell");
	verticesOnCellSrc = new NCField<int>(globalMeshFile, "verticesOnCell");
	cellsOnVertexSrc = new NCField<int>(globalMeshFile, "cellsOnVertex");
	cellsOnEdgeSrc = new NCField<int>(globalMeshFile, "cellsOnEdge");
	edgesOnCellSrc = new NCField<int>(globalMeshFile, "edgesOnCell");
	zgridSrc = new NCField<float>(globalMeshFile, "zgrid");
	stop_timer(0, &secs, &nsecs);
	printf("Time to read mesh fields from %s : %i.%9.9i\n", globalMeshFile, secs, nsecs);


	//
	// Create fields to hold zgrid averaged to layer midpoints and then averaged to edges
	// for both the regional and the global meshes
	//
	zmidSrc = new NCField<float>("zmid", 2, "nCells", zgridSrc->dimSize("nCells"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
	zedgeSrc = new NCField<float>("zedge", 2, "nEdges", latEdgeSrc->dimSize("nEdges"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
	zmidDst = new NCField<float>("zmid", 2, "nCells", zgridDst->dimSize("nCells"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);
	zedgeDst = new NCField<float>("zedge", 2, "nEdges", latEdgeDst->dimSize("nEdges"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);


	//
	// Average global grid zgrid field to layer midpoints and then to edges
	//
	start_timer(0);
	zmidSrcArr = zmidSrc->ptr2D();
	zgridSrcArr = zgridSrc->ptr2D();
	avg_to_midpoint((int)zmidSrc->dimSize("nCells"), (int)zmidSrc->dimSize("nVertLevels")+1, zgridSrcArr, zmidSrcArr);
	avg_cell_to_edge((int)latEdgeSrc->dimSize("nEdges"), (int)zmidSrc->dimSize("nVertLevels"), cellsOnEdgeSrc->ptr2D(), zmidSrc->ptr2D(), zedgeSrc->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to average zgridSrc to edges: %i.%9.9i\n", secs, nsecs);


	//
	// Average regional grid zgrid field to layer midpoints and then to edges
	//
	start_timer(0);
	zmidDstArr = zmidDst->ptr2D();
	zgridDstArr = zgridDst->ptr2D();
	avg_to_midpoint((int)zmidDst->dimSize("nCells"), (int)zmidDst->dimSize("nVertLevels")+1, zgridDstArr, zmidDstArr);
	avg_cell_to_edge((int)latEdgeDst->dimSize("nEdges"), (int)zmidDst->dimSize("nVertLevels"), cellsOnEdgeDst->ptr2D(), zmidDst->ptr2D(), zedgeDst->ptr2D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to average zgridDst to edges: %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields on levels (only used for w right now)
	//
	start_timer(0);
	cellLevelMap = new RemapperCell();
	cellLevelMap->computeWeightsCell(latCellDst->dimSize("nCells"), zgridSrc->dimSize("nVertLevelsP1"), zgridDst->dimSize("nVertLevelsP1"), 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), zgridDst->ptr2D(), bdyMaskCellDst->ptr1D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellLevelMap : %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields on layers
	//
	start_timer(0);
	cellLayerMap = new RemapperCell();
	cellLayerMap->computeWeightsCell(latCellDst->dimSize("nCells"), zmidSrc->dimSize("nVertLevels"), zmidDst->dimSize("nVertLevels"), 3,
                                      nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                      latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                      latCellDst->ptr1D(), lonCellDst->ptr1D(), zmidDst->ptr2D(), bdyMaskCellDst->ptr1D());
	stop_timer(0, &secs, &nsecs);
	printf("Time to create cellLayerMap : %i.%9.9i\n", secs, nsecs);


	//
	// Compute weights for remapping cell-based fields to edges (used for interpolating uReconstruct{Zonal,Meridional} to u)
	//
	if (use_reconstruct_winds) {
		start_timer(0);
		cellToEdgeMap = new RemapperCell();
		cellToEdgeMap->computeWeightsCell(latEdgeDst->dimSize("nEdges"), zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"), 3,
                                                  nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                                  latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latVertexSrc->ptr1D(), lonVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                                  latEdgeDst->ptr1D(), lonEdgeDst->ptr1D(), zedgeDst->ptr2D(), bdyMaskCellDst->ptr1D());
		stop_timer(0, &secs, &nsecs);
		printf("Time to create cellToEdgeMap : %i.%9.9i\n", secs, nsecs);
	}


	//
	// Compute weights for remapping edge-based fields on layers (only used for u and v right now; only needed if
	//    --use-reconstruct-winds not specified)
	//
	if (!use_reconstruct_winds) {
		start_timer(0);
		edgeMap = new RemapperEdge();
		edgeMap->computeWeightsEdge(latCellSrc->dimSize("nCells"), latEdgeDst->dimSize("nEdges"),
                                            zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"),
                                            nEdgesOnCellSrc->ptr1D(), cellsOnCellSrc->ptr2D(), edgesOnCellSrc->ptr2D(),
                                            latCellSrc->ptr1D(), lonCellSrc->ptr1D(), latEdgeSrc->ptr1D(), lonEdgeSrc->ptr1D(), zedgeSrc->ptr2D(),
                                            latCellDst->ptr1D(), lonCellDst->ptr1D(), latEdgeDst->ptr1D(), lonEdgeDst->ptr1D(), zedgeDst->ptr2D(), bdyMaskEdgeDst->ptr1D());
		stop_timer(0, &secs, &nsecs);
		printf("Time to create edgeMap : %i.%9.9i\n", secs, nsecs);
	}

	delete bdyMaskCellDst;
	delete bdyMaskEdgeDst;


	//
	// Time-dependent processing for all global input times
	//
	for (itime=argv_idx; itime<argc; itime++) {
		globalFieldFile = argv[itime];


		//
		// Allocate and read global input fields
		//
		start_timer(0);
		try {
			xtime = new NCField<char>(globalFieldFile, "xtime");
		}
		catch (int e) {
			std::cerr << "Error reading xtime field from " << globalFieldFile << std::endl;
			return 1;
		}
		if (use_reconstruct_winds) {
			try {
				uSrc = new NCField<float>(globalFieldFile, "uReconstructZonal");
			}
			catch (int e) {
				std::cerr << "Error reading uReconstructZonal field from " << globalFieldFile << std::endl;
				return 1;
			}
			vSrc = new NCField<float>(globalFieldFile, "uReconstructMeridional");
		}
		else {
			try {
				uSrc = new NCField<float>(globalFieldFile, "u");
			}
			catch (int e) {
				std::cerr << "Error reading u field from " << globalFieldFile << std::endl;
				return 1;
			}
			vSrc = new NCField<float>("v", 3, "Time", (size_t)1, "nEdges", uSrc->dimSize("nEdges"), "nVertLevels", uSrc->dimSize("nVertLevels"));
		}
		thetaSrc = new NCField<float>(globalFieldFile, "theta");
		rhoSrc = new NCField<float>(globalFieldFile, "rho");
		wSrc = new NCField<float>(globalFieldFile, "w");
		stop_timer(0, &secs, &nsecs);
		printf("Time to read time-dependent fields from %s : %i.%9.9i\n", globalFieldFile, secs, nsecs);


		//
		// Allocate fields for interpolated regional fields
		//
		uDst = new NCField<float>("lbc_u", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		vDst = new NCField<float>("lbc_v", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		thetaDst = new NCField<float>("lbc_theta", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		rhoDst = new NCField<float>("lbc_rho", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
		wDst = new NCField<float>("lbc_w", 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevelsP1", zgridDst->dimSize("nVertLevelsP1"));


		uSrcArr = uSrc->ptr3D();
		vSrcArr = vSrc->ptr3D();
		nEdgesOnEdgeSrcArr = nEdgesOnEdgeSrc->ptr1D();
		edgesOnEdgeSrcArr = edgesOnEdgeSrc->ptr2D();
		weightsOnEdgeSrcArr = weightsOnEdgeSrc->ptr2D();
		angleEdgeSrcArr = angleEdgeSrc->ptr1D();


		//
		// Reconstruct the global v field, and rotate the {u,v} vector field so that
		// u is the zonal wind component and v is the meridional wind component
		//
		if (!use_reconstruct_winds) {
			start_timer(0);
			reconstruct_v(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), nEdgesOnEdgeSrcArr, edgesOnEdgeSrcArr, weightsOnEdgeSrcArr, uSrcArr[0], vSrcArr[0]);
			rotate_winds(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), angleEdgeSrcArr, uSrcArr[0], vSrcArr[0], 1);
			stop_timer(0, &secs, &nsecs);
			printf("Time to reconstruct v and rotate winds : %i.%9.9i\n", secs, nsecs);
		}

		rhoSrcArr = rhoSrc->ptr3D();
		rhoDstArr = rhoDst->ptr3D();


		//
		// Set up name of regional output file as lbc.yyyy-mm-dd_hh.nc
		//
		xtimeArr = xtime->ptr2D();
		snprintf(date, (size_t)14, "%s", xtimeArr[0]);
		snprintf(regionalFieldFile, (size_t)64, "lbc.%s.nc", date);


		//
		// Create output file and define fields in it
		//
		stat = nc_create(regionalFieldFile, NC_64BIT_OFFSET, &ncid);

		stat = xtime->defineInFile(ncid);
		stat = uDst->defineInFile(ncid);
		stat = thetaDst->defineInFile(ncid);
		stat = rhoDst->defineInFile(ncid);
		stat = wDst->defineInFile(ncid);

		//
		// Look for scalars to process (qv, qc, qr, etc.)
		//
		for (int i=0; i<NUM_SCALARS; i++) {
			try {
				qxSrc[i] = new NCField<float>(globalFieldFile, qxNames[i]);
				std::cout << "found " << qxNames[i] << " in " << globalFieldFile << std::endl;
				snprintf(qxLbcName, (size_t)64, "lbc_%s", qxNames[i]);
				qxDst[i] = new NCField<float>(qxLbcName, 3, "Time", (size_t)1, "nCells", zmidDst->dimSize("nCells"), "nVertLevels", zmidDst->dimSize("nVertLevels"));
				qxDst[i]->remapFrom(*qxSrc[i], *cellLayerMap);
				stat = qxDst[i]->defineInFile(ncid);
				delete qxSrc[i];
			}
			catch (int e) {
				std::cout << qxNames[i] << " not found in " << globalFieldFile << std::endl;
				qxDst[i] = new NCField<float>();
			}
		}

		stat = nc_enddef(ncid);

		//
		// Interpolate the zonal and meridional winds, and rotate the wind vector field so that u is the normal component
		//
		start_timer(0);
		uDstArr = uDst->ptr3D();
		vDstArr = vDst->ptr3D();
		angleEdgeDstArr = angleEdgeDst->ptr1D();
		if (!use_reconstruct_winds) {
			uDst->remapFrom(*uSrc, *edgeMap);
			vDst->remapFrom(*vSrc, *edgeMap);
		}
		else {
			uDst->remapFrom(*uSrc, *cellToEdgeMap);
			vDst->remapFrom(*vSrc, *cellToEdgeMap);
		}
		rotate_winds(uDst->dimSize("nEdges"), uDst->dimSize("nVertLevels"), angleEdgeDstArr, uDstArr[0], vDstArr[0], 0);


		//
		// Interpolate scalar fields
		//
		thetaDst->remapFrom(*thetaSrc, *cellLayerMap);
		rhoDst->remapFrom(*rhoSrc, *cellLayerMap);
		wDst->remapFrom(*wSrc, *cellLevelMap);
		stop_timer(0, &secs, &nsecs);
		printf("Time to remap fields : %i.%9.9i\n", secs, nsecs);

		start_timer(0);


		//
		// Write interpolated regional fields to output file
		//
		start_timer(0);
		stat = xtime->writeToFile(ncid);
		stat = uDst->writeToFile(ncid);
		stat = thetaDst->writeToFile(ncid);
		stat = rhoDst->writeToFile(ncid);
		stat = wDst->writeToFile(ncid);
		stop_timer(0, &secs, &nsecs);

		for (int i=0; i<NUM_SCALARS; i++) {
			if (qxDst[i]->valid()) {
				stat = qxDst[i]->writeToFile(ncid);
			}
			delete(qxDst[i]);
		}

		printf("Time to write output fields : %i.%9.9i\n", secs, nsecs);

		stat = nc_close(ncid);


		delete xtime;
		delete uSrc;
		delete vSrc;
		delete thetaSrc;
		delete rhoSrc;
		delete wSrc;

		delete uDst;
		delete vDst;
		delete thetaDst;
		delete rhoDst;
		delete wDst;
	}
	

	delete cellsOnEdgeDst;
	delete latCellDst;
	delete lonCellDst;
	delete latEdgeDst;
	delete lonEdgeDst;
	delete angleEdgeDst;
	delete zgridDst;
	delete zedgeDst;
	delete zmidDst;
	delete latCellSrc;
	delete lonCellSrc;
	delete latEdgeSrc;
	delete lonEdgeSrc;
	delete latVertexSrc;
	delete lonVertexSrc;
	delete nEdgesOnCellSrc;
	delete nEdgesOnEdgeSrc;
	delete weightsOnEdgeSrc;
	delete edgesOnEdgeSrc;
	delete angleEdgeSrc;
	delete cellsOnCellSrc;
	delete verticesOnCellSrc;
	delete cellsOnVertexSrc;
	delete cellsOnEdgeSrc;
	delete edgesOnCellSrc;
	delete zgridSrc;
	delete zedgeSrc;
	delete zmidSrc;
	delete cellLayerMap;
	delete cellLevelMap;
	if (!use_reconstruct_winds) {
		delete edgeMap;
	}
	if (use_reconstruct_winds) {
		delete cellToEdgeMap;
	}

	return 0;
}
