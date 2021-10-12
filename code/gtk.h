/* These are helpers for the graphics library, gtk. */
struct gtkContext;

/* Call to start GTK. */
void setupGTK(int *argc, char ***argv, char *title);

/* Set watchtower selected */
void setGTKWatchTower(int vertex);

/* Set vertex selected */
void setGTKVertex(int vertex);

/* Set edge selected */
void setGTKEdge(int edge);

/* Set bounds of GTK - points are scaled and placed based on these values. */
void setupGTKBounds(double newxmin, double newxmax, double newymin, double newymax);

/* Set background points. */
void setupGTKPoints(int pointCount, void *displayData, char *(*getDataString)(void *, int), 
    void (*freeDataString)(char *), double (*getPointx)(int, void *), 
    double (*getPointy)(int, void *));

/* Set DCEL points. */
void setupGTKPolyPoints(int dcelVerticesCount, void *dcel, 
    double (*getDCELVertexX)(void *, int), double (*getDCELVertexY)(void *, int));

/* Set DCEL edges, hasEdge is checked, hasEdgePair only checked if hasEdge is true. */
void setupGTKPolyLines(int edgeCount, void *dcel, 
    int (*getDCELEdgeVertexStart)(void *, int),
    int (*getDCELEdgeVertexPairStart)(void *, int),
    int (*getDCELEdgeVertexEnd)(void *, int),
    int (*getDCELEdgeVertexPairEnd)(void *, int), int (*DCELhasEdge)(void *, int),
    int (*DCELhasEdgePair)(void *, int));

/* Call to set dijkstra interface functions. */
void setupGTKdijkstra(void *digraph, 
    void *(*dijkstra)(void *graph, int source), 
    void (*constructPath)(void *res, void *graph, int *path, 
        int destination), 
    void (*freeDijkstraRes)(void *res));

void setupGTKAction(void (*doAction)(void *, double, double), void *actionContext);

/* Begin interactive GTK portion. */
void runGTKInteraction();

/* Free GTK allocations. */
void freeGTK();
