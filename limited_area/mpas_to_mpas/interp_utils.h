void get_weights_1d(int nSrcLevels, float *srcLevels, float dstLevel, int *nSrcPts, int *srcPts, float *srcWghts);

// Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
// sphere with given radius.
float sphere_distance(float lat1, float lon1, float lat2, float lon2, float radius);

void convert_lx(float *x, float *y, float *z, float radius, float lat, float lon);

int nearest_cell(float target_x, float target_y, float target_z, int start_cell, int nCells,
                 int *nEdgesOnCell, int **cellsOnCell, float *xCell, float *yCell, float *zCell);
int nearest_vertex(float target_x, float target_y, float target_z, int start_vertex,
                   int vertexDegree, int *nEdgesOnCell, int **verticesOnCell, int **cellsOnVertex,
                   float *xCell, float *yCell, float *zCell,
                   float *xVertex, float *yVertex, float *zVertex);

void mpas_wachspress_coordinates(int nVertices, float vertCoords[][3], float *pointInterp, float *mpas_wachspress_coordinates);
float mpas_triangle_signed_area_sphere(float *a, float *b, float *c, float radius);
float max(float a, float b);
float mpas_arc_length(float ax, float ay, float az, float bx, float by, float bz);
float relative_distance(float ax, float ay, float az, float bx, float by, float bz);
