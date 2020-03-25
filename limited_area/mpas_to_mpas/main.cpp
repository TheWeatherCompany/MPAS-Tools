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
    
    NCField<float> *xCellDst, *yCellDst, *zCellDst;
    NCField<float> *xEdgeDst, *yEdgeDst, *zEdgeDst;
    NCField<float> *angleEdgeDst;
    NCField<int> *cellsOnEdgeDst;
    NCField<float> *zgridDst;
    NCField<float> *zmidDst;
    NCField<float> *zedgeDst;
    NCField<float> *uDst;
    NCField<float> *vDst;
    NCField<float> *thetaDst;
    NCField<float> *rhoDst;
    NCField<float> *wDst;
    NCField<float> *presDst;
    NCField<float> *qxDst[NUM_SCALARS];
    NCField<char> *xtime;
    float ***uDstArr;
    float ***vDstArr;
    float **zgridDstArr;
    float **zmidDstArr;
    float *angleEdgeDstArr;
    
    char qxLbcName[64];
    
    const char *globalMeshFile;
    const char *targetMeshFile;
    const char *globalFieldFile;
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
    NCField<float> *uZonalSrc;
    NCField<float> *uMeridionalSrc;
    float **zgridSrcArr;
    float **zmidSrcArr;
    RemapperCell *cellLayerMap;
    RemapperCell *cellLevelMap;
    RemapperCell *cellToEdgeMap;
    RemapperEdge *edgeMap;
    int secs, nsecs, tsecs, tnsecs;
    int itime;
    int argv_idx;
    int use_reconstruct_winds;
    
    
    if (argc < 4) {
        std::cerr << "\nUsage: " << argv[0] << " [--use-reconstruct-winds] <global_IC_file> <target_IC_file> <global_fields_file>+\n\n";
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
    targetMeshFile = argv[argv_idx++];
    
    
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
        angleEdgeDst = new NCField<float>(ncid, "angleEdge");
        cellsOnEdgeDst = new NCField<int>(ncid, "cellsOnEdge");
        zgridDst = new NCField<float>(ncid, "zgrid");
        
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
            xArr[i] = xArr[i]/sphere_radius;
            yArr[i] = yArr[i]/sphere_radius;
            zArr[i] = zArr[i]/sphere_radius;
        }
        size_t nEdges = static_cast<size_t>(xEdgeDst->dimSize("nEdges"));
        xArr = xEdgeDst->ptr1D();
        yArr = yEdgeDst->ptr1D();
        zArr = zEdgeDst->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nEdges; i++) {
            xArr[i] = xArr[i]/sphere_radius;
            yArr[i] = yArr[i]/sphere_radius;
            zArr[i] = zArr[i]/sphere_radius;
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
        if (!use_reconstruct_winds) {
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
        xArr = xCellSrc->ptr1D();
        yArr = yCellSrc->ptr1D();
        zArr = zCellSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nCells; i++) {
            xArr[i] = xArr[i]/sphere_radius;
            yArr[i] = yArr[i]/sphere_radius;
            zArr[i] = zArr[i]/sphere_radius;
        }
        size_t nEdges = static_cast<size_t>(xEdgeSrc->dimSize("nEdges"));
        xArr = xEdgeSrc->ptr1D();
        yArr = yEdgeSrc->ptr1D();
        zArr = zEdgeSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nEdges; i++) {
            xArr[i] = xArr[i]/sphere_radius;
            yArr[i] = yArr[i]/sphere_radius;
            zArr[i] = zArr[i]/sphere_radius;
        }
        size_t nVertices = static_cast<size_t>(xVertexSrc->dimSize("nVertices"));
        xArr = xVertexSrc->ptr1D();
        yArr = yVertexSrc->ptr1D();
        zArr = zVertexSrc->ptr1D();
#pragma omp parallel for simd
        for (size_t i=0; i<nVertices; i++) {
            xArr[i] = xArr[i]/sphere_radius;
            yArr[i] = yArr[i]/sphere_radius;
            zArr[i] = zArr[i]/sphere_radius;
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
    cellLevelMap = new RemapperCell(xCellDst->dimSize("nCells"), zgridSrc->dimSize("nVertLevelsP1"), zgridDst->dimSize("nVertLevelsP1"), 3);
    stop_timer(0, &secs, &nsecs);
    printf("Time to allocate cellLevelMap : %i.%9.9i\n", secs, nsecs);
    tsecs = secs; tnsecs = nsecs;
    start_timer(0);
    cellLevelMap->computeWeightsCell( nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                     xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                     xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                     xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(), zgridDst->ptr2D());
    stop_timer(0, &secs, &nsecs);
    printf("Time to create cellLevelMap : %i.%9.9i\n", secs, nsecs);
    tsecs += secs; tnsecs += nsecs;
    printf("Total Time to create cellLevelMap : %i.%9.9i\n", tsecs, tnsecs);
    
    
    //
    // Compute weights for remapping cell-based fields on layers
    //
    start_timer(0);
    cellLayerMap = new RemapperCell(xCellDst->dimSize("nCells"), zmidSrc->dimSize("nVertLevels"), zmidDst->dimSize("nVertLevels"), 3);
    stop_timer(0, &secs, &nsecs);
    printf("Time to allocate cellLayerMap : %i.%9.9i\n", secs, nsecs);
    tsecs = secs; tnsecs = nsecs;
    start_timer(0);
    cellLayerMap->computeWeightsCell( nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                     xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                     xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                     xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(), zmidDst->ptr2D());
    stop_timer(0, &secs, &nsecs);
    printf("Time to create cellLayerMap : %i.%9.9i\n", secs, nsecs);
    tsecs += secs; tnsecs += nsecs;
    printf("Total Time to create cellLayerMap : %i.%9.9i\n", tsecs, tnsecs);
    
    
    //
    // Compute weights for remapping cell-based fields to edges (used for interpolating uReconstruct{Zonal,Meridional} to u)
    //
    if (use_reconstruct_winds) {
        start_timer(0);
        cellToEdgeMap = new RemapperCell(xEdgeDst->dimSize("nEdges"), zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"), 3);
        stop_timer(0, &secs, &nsecs);
        printf("Time to allocate cellToEdgeMap : %i.%9.9i\n", secs, nsecs);
        tsecs = secs; tnsecs = nsecs;
        start_timer(0);
        cellToEdgeMap->computeWeightsCell(nEdgesOnCellSrc->ptr1D(), verticesOnCellSrc->ptr2D(), cellsOnVertexSrc->ptr2D(),
                                          xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                          xVertexSrc->ptr1D(), yVertexSrc->ptr1D(), zVertexSrc->ptr1D(), zmidSrc->ptr2D(),
                                          xEdgeDst->ptr1D(), yEdgeDst->ptr1D(), zEdgeDst->ptr1D(), zedgeDst->ptr2D());
        // !! zEdgeDst is the z-coordinate of the edge in Cartesian space,
        //    zedgeDst is the height field of the edges
        stop_timer(0, &secs, &nsecs);
        printf("Time to create cellToEdgeMap : %i.%9.9i\n", secs, nsecs);
        tsecs += secs; tnsecs += nsecs;
        printf("Total Time to create cellToEdgeMap : %i.%9.9i\n", tsecs, tnsecs);
    }
    
    
    //
    // Compute weights for remapping edge-based fields on layers (only used for u and v right now; only needed if
    //    --use-reconstruct-winds not specified)
    //
    if (!use_reconstruct_winds) {
        start_timer(0);
        edgeMap = new RemapperEdge(cellsOnCellSrc->dimSize("maxEdges"), xCellSrc->dimSize("nCells"), xEdgeDst->dimSize("nEdges"),
                                   zedgeSrc->dimSize("nVertLevels"), zedgeDst->dimSize("nVertLevels"));
        stop_timer(0, &secs, &nsecs);
        printf("Time to allocate edgeMap : %i.%9.9i\n", secs, nsecs);
        tsecs = secs; tnsecs = nsecs;
        start_timer(0);
        edgeMap->computeWeightsEdge(nEdgesOnCellSrc->ptr1D(), cellsOnCellSrc->ptr2D(), edgesOnCellSrc->ptr2D(),
                                    xCellSrc->ptr1D(), yCellSrc->ptr1D(), zCellSrc->ptr1D(),
                                    xEdgeSrc->ptr1D(), yEdgeSrc->ptr1D(), zEdgeSrc->ptr1D(), zedgeSrc->ptr2D(),
                                    xCellDst->ptr1D(), yCellDst->ptr1D(), zCellDst->ptr1D(),
                                    zEdgeDst->ptr1D(), yEdgeDst->ptr1D(), zEdgeDst->ptr1D(), zedgeDst->ptr2D());
        stop_timer(0, &secs, &nsecs);
        printf("Time to create edgeMap : %i.%9.9i\n", secs, nsecs);
        tsecs += secs; tnsecs += nsecs;
        printf("Total Time to create edgeMap : %i.%9.9i\n", tsecs, tnsecs);
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
        if (use_reconstruct_winds) {
            try {
                uSrc = new NCField<float>(globalFieldFile, "uReconstructZonal");
            }
            catch (int e) {
                std::cerr << "Error reading uReconstructZonal field from " << globalFieldFile << std::endl;
                return 1;
            }
            vSrc = new NCField<float>(globalFieldFile, "uReconstructMeridional");
            uZonalSrc = new NCField<float>(globalFieldFile, "uReconstructZonal");
            uMeridionalSrc = new NCField<float>(globalFieldFile, "uReconstructMeridional");
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
        stop_timer(0, &secs, &nsecs);
        printf("Time to read time-dependent fields from %s : %i.%9.9i\n", globalFieldFile, secs, nsecs);
        
        
        //
        // Allocate fields for interpolated target mesh fields
        //
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
        
        // Reading scalars
        start_timer(0);
        bool scalars_found[NUM_SCALARS];
        stat = nc_open(globalFieldFile, NC_SHARE, &ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
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
        stat = nc_close(ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        stop_timer(0, &secs, &nsecs);
        printf("Finished reading scalars : %i.%9.9i\n", secs, nsecs);
        
        
        //
        // Set up name of target mesh output file as interpolated.yyyy-mm-dd_hh.nc
        //
        start_timer(0);
        xtimeArr = xtime->ptr2D();
        snprintf(date, (size_t)20, "%s", xtimeArr[0]);
        snprintf(targetFieldFile, (size_t)64, "interpolated.%s.nc", date);
        
        //
        // Create output file and define fields in it
        //
        stat = nc_create(targetFieldFile, NC_64BIT_DATA, &ncid);
        for (int i=0; i<NUM_SCALARS; i++) {
            if (scalars_found[i]) {
                stat = qxDst[i]->defineInFile(ncid);
            }
        }
        stat = xtime->defineInFile(ncid);
        stat = uDst->defineInFile(ncid);
        thetaDst = new NCField<float>("theta", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevels", nVertLevels);
        stat = thetaDst->defineInFile(ncid);
        rhoDst = new NCField<float>("rho", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevels", nVertLevels);
        stat = rhoDst->defineInFile(ncid);
        wDst = new NCField<float>("w", 3, "Time", (size_t)1, "nCells", nCells, "nVertLevelsP1", nVertLevels+1);
        stat = wDst->defineInFile(ncid);
        presDst = new NCField<float>("surface_pressure", 2, "Time", (size_t)1, "nCells", nCells);
        stat = presDst->defineInFile(ncid);
        stat = nc_enddef(ncid);
        
        stat = xtime->writeToFile(ncid);
        delete xtime;
        
        stop_timer(0, &secs, &nsecs);
        printf("Completed file definition : %i.%9.9i\n", secs, nsecs);
        
        //
        // Look for scalars to process (qv, qc, qr, etc.)
        //
        start_timer(0);
        int src_ncid;
        stat = nc_open(globalFieldFile, NC_SHARE, &src_ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        for (int i=0; i<NUM_SCALARS; i++) {
            if (scalars_found[i]) {
                NCField<float> * qxSrc = new NCField<float>(src_ncid, qxNames[i]);
                qxDst[i]->remapFrom(*qxSrc, *cellLayerMap);
                stat = qxDst[i]->writeToFile(ncid);
                delete qxSrc;
                delete qxDst[i];
            }
        }
        stat = nc_close(ncid);
        if (stat != NC_NOERR) {
            throw stat;
        }
        stop_timer(0, &secs, &nsecs);
        printf("Time to remap and write scalars : %i.%9.9i\n", secs, nsecs);
        
        start_timer(0);
        
        //
        // Interpolate the zonal and meridional winds, and rotate the wind vector field so that u is the normal component
        //
        if (!use_reconstruct_winds) {
            uDst->remapFrom(*uSrc, *edgeMap);
            vDst->remapFrom(*vSrc, *edgeMap);
            delete edgeMap;
        }
        else {
            uDst->remapFrom(*uSrc, *cellToEdgeMap);
            vDst->remapFrom(*vSrc, *cellToEdgeMap);
            delete cellToEdgeMap;
        }
        uDstArr = uDst->ptr3D();
        vDstArr = vDst->ptr3D();
        angleEdgeDstArr = angleEdgeDst->ptr1D();
        rotate_winds(uDst->dimSize("nEdges"), uDst->dimSize("nVertLevels"), angleEdgeDstArr, uDstArr[0], vDstArr[0], 0);
        stat = uDst->writeToFile(ncid);
        delete uSrc;
        delete vSrc;
        delete uDst;
        delete vDst;
        if (use_reconstruct_winds) {
            delete uZonalSrc;
            delete uMeridionalSrc;
        }
        
        stop_timer(0, &secs, &nsecs);
        printf("Time to remap and write w : %i.%9.9i\n", secs, nsecs);
        
        
        //
        // Interpolate scalar fields
        //
        
        start_timer(0);
        
        printf("Remapping other scalars\n");
        NCField<float> * thetaSrc = new NCField<float>(globalFieldFile, "theta");
        thetaDst->remapFrom(*thetaSrc, *cellLayerMap);
        printf("Writing theta\n");
        stat = thetaDst->writeToFile(ncid);
        delete thetaSrc;
        delete thetaDst;
        
        NCField<float> * rhoSrc = new NCField<float>(globalFieldFile, "rho");
        rhoDst->remapFrom(*rhoSrc, *cellLayerMap);
        printf("Writing rho\n");
        stat = rhoDst->writeToFile(ncid);
        delete rhoSrc;
        delete rhoDst;
        
        NCField<float> * wSrc = new NCField<float>(globalFieldFile, "w");
        wDst->remapFrom(*wSrc, *cellLevelMap);
        printf("Writing w\n");
        stat = wDst->writeToFile(ncid);
        delete wSrc;
        delete wDst;
        
        NCField<float> * presSrc = new NCField<float>(globalFieldFile, "surface_pressure");
        presDst->remapFrom(*presSrc, *cellLayerMap);
        printf("Writing surface_pressure\n");
        stat = presDst->writeToFile(ncid);
        delete presSrc;
        delete presDst;
        
        stop_timer(0, &secs, &nsecs);
        printf("Time to remap and write other fields : %i.%9.9i\n", secs, nsecs);
        
        stat = nc_close(ncid);
        
    }
    
    
    delete angleEdgeDst;
    delete cellLayerMap;
    delete cellLevelMap;
    
    return 0;
}
