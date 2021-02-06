#include <iostream>
#include <stdio.h>
#include <math.h>
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
    
    NCField<float> *xCellDst, *yCellDst, *zCellDst;
    NCField<int> *landmaskDst;
    NCField<float> *xEdgeDst, *yEdgeDst, *zEdgeDst;
    NCField<float> *angleEdgeDst;
    NCField<int> *cellsOnEdgeDst;
    NCField<float> *zgridDst;
    NCField<float> *zmidDst;
    NCField<float> *zsoilLevsDst;
    NCField<float> *zedgeDst;
    NCField<float> *uDst;
    NCField<float> *u10Dst;
    NCField<float> *vDst;
    NCField<float> *v10Dst;
    NCField<float> *q2Dst;
    NCField<float> *t2mDst;
    NCField<float> *thetaDst;
    NCField<float> *rhoDst;
    NCField<float> *wDst;
    NCField<float> *presDst;
    NCField<float> *tmnDst;
    NCField<float> *skintempDst;
    NCField<float> *snowDst;
    NCField<float> *snowcDst;
    NCField<float> *snowhDst;
    NCField<float> *xiceDst;
    NCField<float> *seaiceDst;
    NCField<float> *dzsDst;
    NCField<float> *smoisDst;
    NCField<float> *tslbDst;
    NCField<float> *qxDst[NUM_SCALARS];
    NCField<char> *xtime;
    float **zgridDstArr;
    float **zmidDstArr;
    
    char qxLbcName[64];
    
    const char *globalMeshFile;
    const char *targetMeshFile;
    const char *globalFieldFile;
    const char *outputFile;
    char targetFieldFile[64];
    char **xtimeArr;
    char date[20];
    
    NCField<float> *xCellSrc, *yCellSrc, *zCellSrc;
    NCField<float> *xEdgeSrc, *yEdgeSrc, *zEdgeSrc;
    NCField<float> *xVertexSrc, *yVertexSrc, *zVertexSrc;
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
    float **zgridSrcArr;
    float **zmidSrcArr;
    RemapperCell *cellLayerMap;
    RemapperCell *cellLevelMap;
    RemapperCell *cellToEdgeMap;
    NCField<int> *landmaskSrc;
    RemapperEdge *edgeMap;
    int secs, nsecs;
    int itime;
    int argv_idx;
    int use_reconstruct_winds;
    // The next three options exclude fields in the interpolation.
    //
    // Exclude Theta, Rho, Q[x]
    int exclude_t_qx_rho;
    // Exclude Winds
    int exclude_winds;
    // Exclude T2, Q2, U10, V10, SoilT, SoilM, 
    // Snow, Ice, TMN, PSFC, TSK
    int exclude_lsm_diag;
    
   
    // 5 Required Arguments 
    if (argc < 5) {
        std::cerr << "\nUsage: " << argv[0] << " [--exclude-t-qx-rho] [--exclude-lsm-diag]"
                                               " [--exclude-winds] [--use-reconstruct-winds]"
                                               " <global_mesh_file> <target_IC_file>"
                                               " <output> <global_fields_file>\n\n";
        return 1;
    }
    
    // Look for --exclude-t-qx-rho option. 
    if (strcmp(argv[1], "--exclude-t-qx-rho") == 0) { 
        std::cout << "Excluding Theta, Rho, and Q[x] ...\n";
        argv_idx = 2;
        exclude_t_qx_rho = 1;
    }
    else {
        std::cout << "Including Theta, Rho, and Q[x] ...\n";
        argv_idx = 1;
        exclude_t_qx_rho = 0;
    }

    // Look for --exclude-lsm-diag option.
    if (strcmp(argv[argv_idx], "--exclude-lsm-diag") == 0) {
        std::cout << "Excluding LSM fields and general diagnostics ...\n";
        argv_idx++;
        exclude_lsm_diag = 1;
    }
    else {
        std::cout << "Including LSM fields and general diagnostics ...\n";
        exclude_lsm_diag = 0;
    }

    // Look for --exclude-winds option.
    if (strcmp(argv[argv_idx], "--exclude-winds") == 0) {
        std::cout << "Excluding wind field(s) ...\n";
        argv_idx++;
        exclude_winds = 1;
    }
    else {
        std::cout << "Including wind field(s) ...\n";
        exclude_winds = 0;
    }
    
    // Look for --use-reconstruct-winds option.
    if (strcmp(argv[argv_idx], "--use-reconstruct-winds") == 0) {
        if (!exclude_winds) {
            std::cout << "Using reconstructed winds at cell centers ...\n";
        }
        argv_idx++;
        use_reconstruct_winds = 1;
    }
    else {
        //
        // Try to catch other uses options besides --use-reconstruct-winds
        //
        if (argv[argv_idx][0] == '-') {
            std::cerr << "Unrecognized option: " << argv[argv_idx] << std::endl;
            return 1;
        }
       
        if (!exclude_winds) { 
            std::cout << "Using normal component of winds at edges ...\n";
        }
        use_reconstruct_winds = 0;
    }

    globalMeshFile = argv[argv_idx++];
    targetMeshFile = argv[argv_idx++];
    outputFile = argv[argv_idx++];
  
    //
    // Read mesh description fields from target IC file
    //
    start_timer(0);
    try {
        int stat, ncid;
        stat = nc_open(targetMeshFile, NC_SHARE, &ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        xCellDst = new NCField<float>(ncid, "xCell");
        yCellDst = new NCField<float>(ncid, "yCell");
        zCellDst = new NCField<float>(ncid, "zCell");
        xEdgeDst = new NCField<float>(ncid, "xEdge");
        yEdgeDst = new NCField<float>(ncid, "yEdge");
        zEdgeDst = new NCField<float>(ncid, "zEdge");
        landmaskDst = new NCField<int>(ncid, "landmask");
        angleEdgeDst = new NCField<float>(ncid, "angleEdge");
        cellsOnEdgeDst = new NCField<int>(ncid, "cellsOnEdge");
        zgridDst = new NCField<float>(ncid, "zgrid");
        zsoilLevsDst = new NCField<float>(ncid, "zs");
        
        float sphere_radius;
        stat = nc_get_att_float (ncid, NC_GLOBAL, "sphere_radius", &sphere_radius);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        stat = nc_close(ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        // normalize
        float *xArr, *yArr, *zArr;
        size_t nCells = static_cast<size_t>(xCellDst->dimSize("nCells"));
        xArr = xCellDst->ptr1D();
        yArr = yCellDst->ptr1D();
        zArr = zCellDst->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nCells; i++) {
            float r = sqrtf(xArr[i]*xArr[i] + yArr[i]*yArr[i] + zArr[i]*zArr[i]);
            xArr[i] = xArr[i]/r;
            yArr[i] = yArr[i]/r;
            zArr[i] = zArr[i]/r;
        }
        size_t nEdges = static_cast<size_t>(xEdgeDst->dimSize("nEdges"));
        xArr = xEdgeDst->ptr1D();
        yArr = yEdgeDst->ptr1D();
        zArr = zEdgeDst->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nEdges; i++) {
            float r = sqrtf(xArr[i]*xArr[i] + yArr[i]*yArr[i] + zArr[i]*zArr[i]);
            xArr[i] = xArr[i]/r;
            yArr[i] = yArr[i]/r;
            zArr[i] = zArr[i]/r;
        }
    }
    catch (int e) {
        std::cerr << "Error reading target fields from " << targetMeshFile << std::endl;
        return 1;
    }
    stop_timer(0, &secs, &nsecs);
    printf("Time to read mesh fields from %s : %i.%9.9i\n", targetMeshFile, secs, nsecs);
    
    
    //
    // Read mesh description fields from global IC file
    //
    size_t nCellSrc = 0;
    start_timer(0);
    try {
        int stat, ncid;
        stat = nc_open(globalMeshFile, NC_SHARE, &ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        xCellSrc = new NCField<float>(ncid, "xCell");
        yCellSrc = new NCField<float>(ncid, "yCell");
        zCellSrc = new NCField<float>(ncid, "zCell");
        xEdgeSrc = new NCField<float>(ncid, "xEdge");
        yEdgeSrc = new NCField<float>(ncid, "yEdge");
        zEdgeSrc = new NCField<float>(ncid, "zEdge");
        xVertexSrc = new NCField<float>(ncid, "xVertex");
        yVertexSrc = new NCField<float>(ncid, "yVertex");
        zVertexSrc = new NCField<float>(ncid, "zVertex");
        nEdgesOnCellSrc = new NCField<int>(ncid, "nEdgesOnCell");
        if (!exclude_winds and !use_reconstruct_winds) {
            nEdgesOnEdgeSrc = new NCField<int>(ncid, "nEdgesOnEdge");
            weightsOnEdgeSrc = new NCField<float>(ncid, "weightsOnEdge");
            edgesOnEdgeSrc = new NCField<int>(ncid, "edgesOnEdge");
            angleEdgeSrc = new NCField<float>(ncid, "angleEdge");
        }
        cellsOnCellSrc = new NCField<int>(ncid, "cellsOnCell");
        verticesOnCellSrc = new NCField<int>(ncid, "verticesOnCell");
        cellsOnVertexSrc = new NCField<int>(ncid, "cellsOnVertex");
        cellsOnEdgeSrc = new NCField<int>(ncid, "cellsOnEdge");
        edgesOnCellSrc = new NCField<int>(ncid, "edgesOnCell");
        zgridSrc = new NCField<float>(ncid, "zgrid");
        landmaskSrc = new NCField<int>(ncid, "landmask");
 
        float sphere_radius;
        stat = nc_get_att_float (ncid, NC_GLOBAL, "sphere_radius", &sphere_radius);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        stat = nc_close(ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        // normalize
        float *xArr, *yArr, *zArr;
        size_t nCells = static_cast<size_t>(xCellSrc->dimSize("nCells"));
        nCellSrc = nCells;
        fprintf(stderr,"nCells from xCellSrc: %d\n",nCells);
        xArr = xCellSrc->ptr1D();
        yArr = yCellSrc->ptr1D();
        zArr = zCellSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nCells; i++) {
            float r = sqrtf(xArr[i]*xArr[i] + yArr[i]*yArr[i] + zArr[i]*zArr[i]);
            xArr[i] = xArr[i]/r;
            yArr[i] = yArr[i]/r;
            zArr[i] = zArr[i]/r;
        }
        size_t nEdges = static_cast<size_t>(xEdgeSrc->dimSize("nEdges"));
        xArr = xEdgeSrc->ptr1D();
        yArr = yEdgeSrc->ptr1D();
        zArr = zEdgeSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nEdges; i++) {
            float r = sqrtf(xArr[i]*xArr[i] + yArr[i]*yArr[i] + zArr[i]*zArr[i]);
            xArr[i] = xArr[i]/r;
            yArr[i] = yArr[i]/r;
            zArr[i] = zArr[i]/r;
        }
        size_t nVertices = static_cast<size_t>(xVertexSrc->dimSize("nVertices"));
        xArr = xVertexSrc->ptr1D();
        yArr = yVertexSrc->ptr1D();
        zArr = zVertexSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nVertices; i++) {
            float r = sqrtf(xArr[i]*xArr[i] + yArr[i]*yArr[i] + zArr[i]*zArr[i]);
            xArr[i] = xArr[i]/r;
            yArr[i] = yArr[i]/r;
            zArr[i] = zArr[i]/r;
        }
    }
    catch (int e) {
        std::cerr << "Error reading global mesh fields from " << globalMeshFile << std::endl;
        return 1;
    }
    stop_timer(0, &secs, &nsecs);
    printf("Time to read mesh fields from %s : %i.%9.9i\n", globalMeshFile, secs, nsecs);
    
    
    //
    // Create fields to hold zgrid averaged to layer midpoints and then averaged to edges
    // for both the target global and the input global meshes
    //
    zmidSrc = new NCField<float>("zmid", 2, "nCells", zgridSrc->dimSize("nCells"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
    zedgeSrc = new NCField<float>("zedge", 2, "nEdges", xEdgeSrc->dimSize("nEdges"), "nVertLevels", zgridSrc->dimSize("nVertLevelsP1")-1);
    zmidDst = new NCField<float>("zmid", 2, "nCells", zgridDst->dimSize("nCells"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);
    zedgeDst = new NCField<float>("zedge", 2, "nEdges", xEdgeDst->dimSize("nEdges"), "nVertLevels", zgridDst->dimSize("nVertLevelsP1")-1);
    
    //
    // Average global grid zgrid field to layer midpoints and then to edges
    //
    start_timer(0);
    zmidSrcArr = zmidSrc->ptr2D();
    zgridSrcArr = zgridSrc->ptr2D();
    avg_to_midpoint((int)zmidSrc->dimSize("nCells"), (int)zmidSrc->dimSize("nVertLevels")+1, zgridSrcArr, zmidSrcArr);
    avg_cell_to_edge((int)xEdgeSrc->dimSize("nEdges"), (int)zmidSrc->dimSize("nVertLevels"), cellsOnEdgeSrc->ptr2D(), zmidSrc->ptr2D(), zedgeSrc->ptr2D());
    stop_timer(0, &secs, &nsecs);
    printf("Time to average zgridSrc to edges: %i.%9.9i\n", secs, nsecs);
    
    
    //
    // Average target grid zgrid field to layer midpoints and then to edges
    //
    start_timer(0);
    zmidDstArr = zmidDst->ptr2D();
    zgridDstArr = zgridDst->ptr2D();
    avg_to_midpoint((int)zmidDst->dimSize("nCells"), (int)zmidDst->dimSize("nVertLevels")+1, zgridDstArr, zmidDstArr);
    avg_cell_to_edge((int)xEdgeDst->dimSize("nEdges"), (int)zmidDst->dimSize("nVertLevels"), cellsOnEdgeDst->ptr2D(), zmidDst->ptr2D(), zedgeDst->ptr2D());
    stop_timer(0, &secs, &nsecs);
    printf("Time to average zgridDst to edges: %i.%9.9i\n", secs, nsecs);
    
    
    //
    // Compute weights for remapping cell-based fields on levels (only used for w right now)
    //
    start_timer(0);
    cellLevelMap = new RemapperCell(xCellDst->dimSize("nCells"), zgridSrc->dimSize("nVertLevelsP1"), zgridDst->dimSize("nVertLevelsP1"), 3, zsoilLevsDst->dimSize("nSoilLevels"));
    cellLevelMap->computeWeightsCell( nCellSrc, nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(), cellsOnCellSrc->ptr2D(),
                                     xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                     xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                     landmaskSrc->ptr1D(), 
                                     xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(), landmaskDst->ptr1D(), 
                                     zgridDst->ptr2D() );
    stop_timer(0, &secs, &nsecs);
    printf("Time to create cellLevelMap : %i.%9.9i\n", secs, nsecs);
    
    

    //
    // Compute weights for remapping cell-based fields on layers
    //
    start_timer(0);
    cellLayerMap = new RemapperCell(xCellDst->dimSize("nCells"), zmidSrc->dimSize("nVertLevels"), zmidDst->dimSize("nVertLevels"), 3, zsoilLevsDst->dimSize("nSoilLevels"));
    cellLayerMap->computeWeightsCell( nCellSrc, nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(), cellsOnCellSrc->ptr2D(),
                                     xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                     xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                     landmaskSrc->ptr1D(),
                                     xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(), landmaskDst->ptr1D(),
                                     zmidDst->ptr2D() );
    stop_timer(0, &secs, &nsecs);
    printf("Time to create cellLayerMap : %i.%9.9i\n", secs, nsecs);

    //
    // Compute weights for remapping cell-based fields to edges (used for interpolating uReconstruct{Zonal,Meridional} to u)
    //
    if (!exclude_winds and use_reconstruct_winds) {
        start_timer(0);
        cellToEdgeMap = new RemapperCell(xEdgeDst->dimSize("nEdges"), zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"), 3);
        cellToEdgeMap->computeWeightsCell(nCellSrc, nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(), cellsOnCellSrc->ptr2D(),
                                          xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                          xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                          landmaskSrc->ptr1D(),
                                          xEdgeDst->ptr1D(), yEdgeDst->ptr1D(), zEdgeDst->ptr1D(), landmaskDst->ptr1D(),
                                          zedgeDst->ptr2D() );
        // !! zEdgeDst is the z-coordinate of the edge in Cartesian space,
        //    zedgeDst is the height field of the edges
        stop_timer(0, &secs, &nsecs);
        printf("Time to create cellToEdgeMap : %i.%9.9i\n", secs, nsecs);
    }
    
    
    //
    // Compute weights for remapping edge-based fields on layers (only used for u and v right now; only needed if
    //    --use-reconstruct-winds not specified)
    //
    if (!exclude_winds and !use_reconstruct_winds) {
        start_timer(0);
        edgeMap = new RemapperEdge(cellsOnCellSrc->dimSize("maxEdges"), xCellSrc->dimSize("nCells"), xEdgeDst->dimSize("nEdges"),
                                   zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"));
        edgeMap->computeWeightsEdge(nEdgesOnCellSrc->ptr1D(), cellsOnCellSrc->ptr2D(), edgesOnCellSrc->ptr2D(),
                                    xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                    xEdgeSrc->ptr1D(), yEdgeSrc->ptr1D(), zEdgeSrc->ptr1D(), zedgeSrc->ptr2D(),
                                    xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(),
                                    zEdgeDst->ptr1D(), yEdgeDst->ptr1D(), zEdgeDst->ptr1D(), zedgeDst->ptr2D());
        stop_timer(0, &secs, &nsecs);
        printf("Time to create edgeMap : %i.%9.9i\n", secs, nsecs);
    }
    
    delete xCellDst, yCellDst, zCellDst;
    delete xEdgeDst, yEdgeDst, zEdgeDst;
    delete xCellSrc, yCellSrc, zCellSrc;
    delete xEdgeSrc, yEdgeSrc, zEdgeSrc;
    delete xVertexSrc, yVertexSrc, zVertexSrc;
    delete zgridDst;
    delete cellsOnEdgeDst, edgesOnCellSrc;
    delete cellsOnCellSrc;
    delete verticesOnCellSrc, cellsOnVertexSrc;
    delete cellsOnEdgeSrc;
    delete nEdgesOnCellSrc;
    delete zmidSrc;
    delete zgridSrc;
    delete zedgeSrc, zedgeDst;
    
    size_t nCells = zmidDst->dimSize("nCells");
    size_t nVertLevels = zmidDst->dimSize("nVertLevels");
    size_t nSoilLevels = zsoilLevsDst->dimSize("nSoilLevels");
    fprintf(stderr,"nSoilLevels: %d nVertLevels: %d\n", nSoilLevels,nVertLevels);
    delete zsoilLevsDst; 
    delete zmidDst;
    
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
        if (!exclude_winds) {
            if (use_reconstruct_winds) {
                // Defer reading uReconstructZonal/Meridional till needed
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
        }
        stop_timer(0, &secs, &nsecs);
        printf("Time to read time-dependent fields from %s : %i.%9.9i\n", globalFieldFile, secs, nsecs);
        
        
        //
        // Allocate fields for interpolated target mesh fields
        //
        if (!exclude_winds) {

            uDst = new NCField<float>("u", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", nVertLevels);
            vDst = new NCField<float>("v", 3, "Time", (size_t)1, "nEdges", angleEdgeDst->dimSize("nEdges"), "nVertLevels", nVertLevels);
        
            //
            // Reconstruct the global v field, and rotate the {u,v} vector field so that
            // u is the zonal wind component and v is the meridional wind component
            //
            if (!use_reconstruct_winds) {
                float ***uSrcArr = uSrc->ptr3D();
                float ***vSrcArr = vSrc->ptr3D();
                int * nEdgesOnEdgeSrcArr = nEdgesOnEdgeSrc->ptr1D();
                int ** edgesOnEdgeSrcArr = edgesOnEdgeSrc->ptr2D();
                float ** weightsOnEdgeSrcArr = weightsOnEdgeSrc->ptr2D();
                float * angleEdgeSrcArr = angleEdgeSrc->ptr1D();
                start_timer(0);
                reconstruct_v(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), nEdgesOnEdgeSrcArr, edgesOnEdgeSrcArr, weightsOnEdgeSrcArr, uSrcArr[0], vSrcArr[0]);
                rotate_winds(uSrc->dimSize("nEdges"), uSrc->dimSize("nVertLevels"), angleEdgeSrcArr, uSrcArr[0], vSrcArr[0], 1);
                stop_timer(0, &secs, &nsecs);
                printf("Time to reconstruct v and rotate winds : %i.%9.9i\n", secs, nsecs);
                delete nEdgesOnEdgeSrc;
                delete weightsOnEdgeSrc;
                delete edgesOnEdgeSrc;
                delete angleEdgeSrc;
            }
    
        }
    
        // Identifying which scalar qx's are available from src
        bool scalars_found[NUM_SCALARS];
        stat = nc_open(globalFieldFile, NC_SHARE, &ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }

        if (!exclude_t_qx_rho) {
            for (int i=0; i<NUM_SCALARS; i++) {
                try {
                    int varid;
                    stat = nc_inq_varid(ncid, qxNames[i], &varid);
                    if (stat != NC_NOERR) throw stat;
                    std::cout << "found " << qxNames[i] << " in " << globalFieldFile << std::endl;
                    snprintf(qxLbcName, (size_t)64, "%s", qxNames[i]);
                    qxDst[i] = new NCField<float>(qxLbcName, 3, "Time", (size_t)1, "nCells", nCells, "nVertLevels", nVertLevels);
                    scalars_found[i] = true;
                } catch (int e) {
                    std::cout << qxNames[i] << " not found in " << globalFieldFile << std::endl;
                    qxDst[i] = new NCField<float>();
                    scalars_found[i] = false;
                }
            }
        }

        stat = nc_close(ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        
        //
        // Set up name of target mesh output file as interpolated.yyyy-mm-dd_hh.nc
        //
        start_timer(0);
        xtimeArr = xtime->ptr2D();
        snprintf(date, (size_t)20, "%s", xtimeArr[0]);
        snprintf(targetFieldFile, (size_t)64, outputFile, date);
        
        //
        // Create output file and define fields in it
        //
        stat = nc_create(targetFieldFile, NC_64BIT_DATA|NC_CLOBBER, &ncid);
        stat = xtime->defineInFile(ncid);
        if (!exclude_t_qx_rho) {
            for (int i=0; i<NUM_SCALARS; i++) {
                if (scalars_found[i]) {
                    stat = qxDst[i]->defineInFile(ncid);
                }
            }
        }

        if (!exclude_winds) {
            stat = uDst->defineInFile(ncid);
        } 

        if (!exclude_t_qx_rho) {
            thetaDst = new NCField<float>("theta", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevels", nVertLevels);
            stat = thetaDst->defineInFile(ncid);
            rhoDst = new NCField<float>("rho", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevels", nVertLevels);
            stat = rhoDst->defineInFile(ncid);
        }

        if (!exclude_winds) {
            wDst = new NCField<float>("w", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevelsP1", nVertLevels+1);
            stat = wDst->defineInFile(ncid);
        }

        // Added by Todd Hutchinson during upgrade to be able to use ecmwf input
        //
        if (!exclude_lsm_diag) {
        
            tmnDst = new NCField<float>("tmn", 2, "Time", (size_t)1, "nCells", nCells);
            stat = tmnDst->defineInFile(ncid);

            skintempDst = new NCField<float>("skintemp", 2, "Time", (size_t)1, "nCells", nCells);
            stat = skintempDst->defineInFile(ncid);

            snowDst = new NCField<float>("snow", 2, "Time", (size_t)1, "nCells", nCells);
            stat = snowDst->defineInFile(ncid);

            snowcDst = new NCField<float>("snowc", 2, "Time", (size_t)1, "nCells", nCells);
            stat = snowcDst->defineInFile(ncid);

            snowhDst = new NCField<float>("snowh", 2, "Time", (size_t)1, "nCells", nCells);
            stat = snowhDst->defineInFile(ncid);

            xiceDst = new NCField<float>("xice", 2, "Time", (size_t)1, "nCells", nCells);
            stat = xiceDst->defineInFile(ncid);

            seaiceDst = new NCField<float>("seaice", 2, "Time", (size_t)1, "nCells", nCells);
            stat = seaiceDst->defineInFile(ncid);

            smoisDst = new NCField<float>("smois", 3, "Time", (size_t)1, "nCells", nCells, "nSoilLevels", nSoilLevels);
            stat = smoisDst->defineInFile(ncid);

            tslbDst = new NCField<float>("tslb", 3, "Time", (size_t)1, "nCells", nCells, "nSoilLevels", nSoilLevels);
            stat = tslbDst->defineInFile(ncid);

            dzsDst = new NCField<float>("dzs", 3, "Time", (size_t)1, "nCells", nCells, "nSoilLevels", nSoilLevels);
            stat = dzsDst->defineInFile(ncid);

            u10Dst = new NCField<float>("u10", 2, "Time", (size_t)1, "nCells", nCells);
            stat = u10Dst->defineInFile(ncid);

            v10Dst = new NCField<float>("v10", 2, "Time", (size_t)1, "nCells", nCells);
            stat = v10Dst->defineInFile(ncid);

            q2Dst = new NCField<float>("q2", 2, "Time", (size_t)1, "nCells", nCells);
            stat = q2Dst->defineInFile(ncid);

            t2mDst = new NCField<float>("t2m", 2, "Time", (size_t)1, "nCells", nCells);
            stat = t2mDst->defineInFile(ncid);

            presDst = new NCField<float>("surface_pressure", 2, "Time", (size_t)1, "nCells", nCells);
            stat = presDst->defineInFile(ncid);

        }    

        stat = nc_enddef(ncid);
        
        stat = xtime->writeToFile(ncid);
        if (stat != NC_NOERR) throw stat;
        delete xtime;
        
        stop_timer(0, &secs, &nsecs);
        printf("Completed file definition : %i.%9.9i\n", secs, nsecs);
        
        //
        // Look for scalars to process (qv, qc, qr, etc.)
        //
        int src_ncid;
        stat = nc_open(globalFieldFile, NC_SHARE, &src_ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        if (!exclude_t_qx_rho) {
            for (int i=0; i<NUM_SCALARS; i++) {
                if (scalars_found[i]) {
                    NCField<float> * qxSrc = new NCField<float>(src_ncid, qxNames[i]);
                    start_timer(0);
                    qxDst[i]->remapFrom(*qxSrc, *cellLayerMap, RemapperCell::barycentric);
                    stop_timer(0, &secs, &nsecs);
                    printf("Time to remap qxName : %s %i.%9.9i\n", qxNames[i], secs, nsecs);
                    start_timer(0);
                    stat = qxDst[i]->writeToFile(ncid);
                    stop_timer(0, &secs, &nsecs);
                    printf("Time to write qxName : %s %i.%9.9i\n", qxNames[i], secs, nsecs);
                    if (stat != NC_NOERR) throw stat;
                    delete qxSrc;
                    delete qxDst[i];
                }
            }
        }
        stat = nc_close(src_ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        
        start_timer(0);
        //
        // Interpolate the zonal and meridional winds, and rotate the wind vector field so that u is the normal component
        //
        if (!exclude_winds) {
            start_timer(0);
            if (!use_reconstruct_winds) {
                uDst->remapFrom(*uSrc, *edgeMap);
                vDst->remapFrom(*vSrc, *edgeMap);
                delete uSrc;
                delete vSrc;
            } else {
                uSrc = new NCField<float>(globalFieldFile, "uReconstructZonal");
                uDst->remapFrom(*uSrc, *cellToEdgeMap);
                delete uSrc;
                vSrc = new NCField<float>(globalFieldFile, "uReconstructMeridional");
                vDst->remapFrom(*vSrc, *cellToEdgeMap);
                delete vSrc;
            }
            float *** uDstArr = uDst->ptr3D();
            float *** vDstArr = vDst->ptr3D();
            float * angleEdgeDstArr = angleEdgeDst->ptr1D();
            rotate_winds(uDst->dimSize("nEdges"), uDst->dimSize("nVertLevels"), angleEdgeDstArr, uDstArr[0], vDstArr[0], 0);
            stat = uDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete uDst;
            delete vDst;
        
            stop_timer(0, &secs, &nsecs);
            printf("Time to remap and write winds : %i.%9.9i\n", secs, nsecs);
        }            

        
        //
        // Interpolate scalar fields
        //
        
        if (!exclude_t_qx_rho) {

            printf("Remapping theta, rho, and w\n");

            start_timer(0);
            NCField<float> * thetaSrc = new NCField<float>(globalFieldFile, "theta");
            thetaDst->remapFrom(*thetaSrc, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing theta\n");
            stat = thetaDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete thetaSrc;
            delete thetaDst;
            stop_timer(0, &secs, &nsecs);
            printf("Time to remap and write theta : %i.%9.9i\n", secs, nsecs);

            start_timer(0);
            NCField<float> * rhoSrc = new NCField<float>(globalFieldFile, "rho");
            rhoDst->remapFrom(*rhoSrc, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing rho\n");
            stat = rhoDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete rhoSrc;
            delete rhoDst;
            stop_timer(0, &secs, &nsecs);
            printf("Time to remap and write rho : %i.%9.9i\n", secs, nsecs);
        
        }

        if (!exclude_winds) {
            start_timer(0);
            NCField<float> * wSrc = new NCField<float>(globalFieldFile, "w");
            wDst->remapFrom(*wSrc, *cellLevelMap, RemapperCell::barycentric);
            printf("Writing w\n");
            stat = wDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete wSrc;
            delete wDst;
            stop_timer(0, &secs, &nsecs);
            printf("Time to remap and write w : %i.%9.9i\n", secs, nsecs);
        }

        if (!exclude_lsm_diag) {

            printf("Remapping tmn, skintemp, snow, snowc, snowh, xice, seaice, smois, tslb, etc.\n");
            start_timer(0);

            NCField<float> * tmnSrc = new NCField<float>(globalFieldFile, "tmn");
            tmnDst->remapFrom(*tmnSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing tmn\n");
            stat = tmnDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete tmnSrc;
            delete tmnDst;

            NCField<float> * skintempSrc = new NCField<float>(globalFieldFile, "skintemp");
            skintempDst->remapFrom(*skintempSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing skintemp\n");
            stat = skintempDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete skintempSrc;
            delete skintempDst;

            NCField<float> * snowSrc = new NCField<float>(globalFieldFile, "snow");
            snowDst->remapFrom(*snowSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing snow\n");
            stat = snowDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete snowSrc;
            delete snowDst;

            NCField<float> * snowcSrc = new NCField<float>(globalFieldFile, "snowc");
            snowcDst->remapFrom(*snowcSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing snowc\n");
            stat = snowcDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete snowcSrc;
            delete snowcDst;

            NCField<float> * snowhSrc = new NCField<float>(globalFieldFile, "snowh");
            snowhDst->remapFrom(*snowhSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing snowh\n");
            stat = snowhDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete snowhSrc;
            delete snowhDst;

            NCField<float> * xiceSrc = new NCField<float>(globalFieldFile, "xice");
            xiceDst->remapFrom(*xiceSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing xice\n");
            stat = xiceDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete xiceSrc;
            delete xiceDst;
        
            NCField<float> * seaiceSrc = new NCField<float>(globalFieldFile, "seaice");
            seaiceDst->remapFrom(*seaiceSrc, *cellLayerMap, RemapperCell::nearest_samelandmask);
            printf("Writing seaice\n");
            stat = seaiceDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete seaiceSrc;
            delete seaiceDst;
        
            NCField<float> * smoisSrc = new NCField<float>(globalFieldFile, "smois");
            smoisDst->remapFrom(*smoisSrc, *cellLayerMap, RemapperCell::nearest_samelandmask_soil);
            printf("Writing smois\n");
            stat = smoisDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete smoisSrc;
            delete smoisDst;

            NCField<float> * tslbSrc = new NCField<float>(globalFieldFile, "tslb");
            tslbDst->remapFrom(*tslbSrc, *cellLayerMap, RemapperCell::nearest_samelandmask_soil);
            printf("Writing tslb\n");
            stat = tslbDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete tslbSrc;
            delete tslbDst;

            NCField<float> * dzsSrc = new NCField<float>(globalFieldFile, "dzs");
            dzsDst->remapFrom(*dzsSrc, *cellLayerMap, RemapperCell::nearest_samelandmask_soil);
            printf("Writing dzs\n");
            stat = dzsDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete dzsSrc;
            delete dzsDst;

            NCField<float> * u10Src = new NCField<float>(globalFieldFile, "u10");
            u10Dst->remapFrom(*u10Src, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing u10\n");
            stat = u10Dst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete u10Src;
            delete u10Dst;

            NCField<float> * v10Src = new NCField<float>(globalFieldFile, "v10");
            v10Dst->remapFrom(*v10Src, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing v10\n");
            stat = v10Dst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete v10Src;
            delete v10Dst;

            NCField<float> * q2Src = new NCField<float>(globalFieldFile, "q2");
            q2Dst->remapFrom(*q2Src, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing q2\n");
            stat = q2Dst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete q2Src;
            delete q2Dst;

            NCField<float> * t2mSrc = new NCField<float>(globalFieldFile, "t2m");
            t2mDst->remapFrom(*t2mSrc, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing t2m\n");
            stat = t2mDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete t2mSrc;
            delete t2mDst;

            NCField<float> * presSrc = new NCField<float>(globalFieldFile, "surface_pressure");
            presDst->remapFrom(*presSrc, *cellLayerMap, RemapperCell::barycentric);
            printf("Writing surface_pressure\n");
            stat = presDst->writeToFile(ncid);
            if (stat != NC_NOERR) throw stat;
            delete presSrc;
            delete presDst;
        
            stop_timer(0, &secs, &nsecs);
            printf("Time to remap and write LSM fields and general diagnostics : %i.%9.9i\n", secs, nsecs);
 
        }
        
        stat = nc_close(ncid);
        
    }
    
    if (!exclude_winds) {
        if (use_reconstruct_winds) {
            delete cellToEdgeMap;
        }
        else {
            delete edgeMap;
        }
    }
        
    delete angleEdgeDst;
    delete cellLayerMap;
    delete cellLevelMap;
    
    return 0;
}
