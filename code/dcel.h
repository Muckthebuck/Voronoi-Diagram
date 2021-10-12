#include <stdio.h>
#include "watchtowerStruct.h"
struct face;
struct split;
struct bisector;
struct intersection;

/* Reads the polygon from the given file. */
struct DCEL *readPolygonFile(char *polygonfileName);

/* Reads the next split from the given file. */
struct split *readNextSplit(FILE *splitfile);

/* Frees a given split.  */
void freeSplit(struct split *split);

/* Reads the next bisector from the given file. */
struct bisector *readNextBisector(FILE *bisectorfile);

/* Gets the string for the given bisector equation. */
char *getBisectorEquation(struct bisector *b);

/* Frees the given bisector. */
void freeBisector(struct bisector *bisector);

/* Representation of no face */
#define NOFACE (-1)

/* Default face for intersections. */
#define DEFAULT_FACE 0
/* Default minimum length for bisector in each direction */
#define DEFAULTMINLENGTH (200)
/* 
    Gets the intersection between the given bisector and the given DCEL for the given 
    face. 
*/
struct intersection *getIntersection(struct bisector *b, struct DCEL *dcel, int face,
    double minLength);
/* Gets the string for the given intersection. */
char *getIntersectionString(struct intersection *intersection);
/* Frees a given intersection. */
void freeIntersection(struct intersection *intersection);

/* Applies a given split to the DCEL. */
void applySplit(struct split *split, struct DCEL *dcel);

/* Frees the given DCEL */
void freeDCEL(struct DCEL *dcel);

/* Gets the number of faces in the DCEL. */
int getFaceCount(struct DCEL *dcel);

/* Returns 1 if the given x,y point is inside the given face. */
int inFace(struct DCEL *dcel, double x, double y, int faceIndex);

/* Gets the diameter of the given face. */
double getDiameter(struct DCEL *dcel, int faceIndex);

/*
    Adds the watchtower to the Voronoi diagram represented by the given DCEL,
    applying required splits and setting the watchtower as required.
*/
void incrementalVoronoi(struct DCEL *dcel, struct watchtowerStruct *watchTower);

/* Returns the number of vertices in the DCEL. */
int getDCELPointCount(struct DCEL *dcel);

/* Get x value of given vertex in DCEL. */
double getDCELVertexX(struct DCEL *dcel, int vertex);

/* Get y value of given vertex in DCEL. */
double getDCELVertexY(struct DCEL *dcel, int vertex);

/* Returns the number of edges in the DCEL. */
int getDCELEdgeCount(struct DCEL *dcel);

/* Get start vertex of given edge in DCEL. */
int getDCELEdgeVertexStart(struct DCEL *dcel, int edge);

/* Get end vertex of given edge in DCEL. */
int getDCELEdgeVertexEnd(struct DCEL *dcel, int edge);

/* Get start vertex of paired given edge in DCEL. */
int getDCELEdgeVertexPairStart(struct DCEL *dcel, int edge);

/* Get end vertex of paired given edge in DCEL. */
int getDCELEdgeVertexPairEnd(struct DCEL *dcel, int edge);

/* Check if the DCEL has the given edge. */
int DCELhasEdge(struct DCEL *dcel, int edge);

/* Check if the DCEL has a pair for the given edge. */
int DCELhasEdgePair(struct DCEL *dcel, int edge);
