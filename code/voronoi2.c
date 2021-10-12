/*
    Written by Grady Fitzpatrick for COMP20003 Assignment 2
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "readData.h"
#include "watchtowerStruct.h"
#include "dcel.h"
#include "output.h"
#include "problem.h"

#ifdef USE_GUI
/* GUI Application Extensions */
#include "gtk.h"
#endif

#ifdef USE_GUI
/* Mouse-based action. */
void doAction(void *data, double mouseX, double mouseY);

void doAction(void *data, double mouseX, double mouseY){
    printf("Mouse clicked at (%f, %f)\n", mouseX, mouseY);
}

#define WT_CHANGE 0
#define POLY_CHANGE 1
#define EDGE_CHANGE 2
void doCyleOnClick(int *data, double mouseX, double mouseY);

void doCyleOnClick(int *data, double mouseX, double mouseY){
    data[1] = data[1] + 1;
    if(data[1] % 2 == 1){
        switch(data[0]){
            case WT_CHANGE:
                setGTKWatchTower(data[1] / 2);
                printf("Switched WatchTower to %d\n", data[1] / 2);
                break;
            
            case POLY_CHANGE:
                setGTKVertex(data[1] / 2);
                printf("Switched Vertex to %d\n", data[1] / 2);
                break;
            
            case EDGE_CHANGE:
                setGTKEdge(data[1] / 2);
                printf("Switched Edge to %d\n", data[1] / 2);
                break;
            default:
                break;
        }
    }

    //printf("Mouse clicked at (%f, %f)\n", mouseX, mouseY);
}

#endif

int main(int argc, char **argv){
    char *datafile;
    char *polygonfile;
    char *outputfile;
    char *pointpairfile;
    enum assignmentStage stage = UNSET;
    struct watchtowerStruct **watchTowers;
    int watchTowerCount;
    struct DCEL *dcel;
    struct split *currentSplit = NULL;
    struct bisector *currentBisector = NULL;
    FILE *ppFile = NULL;
    FILE *outFile = NULL;
    
    #ifdef USE_GUI
    /* Use GUI Application Extensions */
    setupGTK(&argc, &argv, "Assignment 2 - Graph Visualisation");
    
    /* Visualisation: Set Bounds */
    #define VICTORIAXMIN (140.9)
    #define VICTORIAXMAX (150.0)
    #define VICTORIAYMIN (-41.85)
    // #define VICTORIAYMIN (-39.2)
    #define VICTORIAYMAX (-32.0)
    //#define VICTORIAYMAX (-33.9)
    
    setupGTKBounds(VICTORIAXMIN, VICTORIAXMAX, VICTORIAYMIN, VICTORIAYMAX);
    #endif
    
    if(argc < 2){
        printArgError("Incorrect argument count", argc, stage);
    }
    
    stage = getArgumentMode(argv[1]);
    if(stage == STAGE_ERROR){
        fprintf(stderr, "Invalid Stage (%s)\n\n", argv[1]);
        printArgError("Stage was invalid\n", argc, stage);
    }
    if(! checkArgCount(stage, argc)){
        printArgError("Incorrect arguments for stage", argc, stage);
    }
    
    switch(stage){
        case STAGE_0:
            datafile = argv[2];
            polygonfile = argv[3];
            outputfile = argv[4];
            watchTowers = readDataFile(datafile, &watchTowerCount);
            dcel = readPolygonFile(polygonfile);
            while((currentSplit = readNextSplit(stdin)) != NULL){
                applySplit(currentSplit, dcel);
                freeSplit(currentSplit);
            }
            #ifdef USE_GUI
            /* Visualisation: Points */
            int pointCount = watchTowerCount;
            setupGTKPoints(pointCount, 
                (void *) watchTowers,
                (char *(*)(void *, int)) &getWTDataString, 
                (void (*)(char *)) &freeWTDataString,
                (double (*)(int, void *)) &getWTx,
                (double (*)(int, void *)) &getWTy
            );

            /* Visualisation: DCEL points */
            int dcelVerticesCount = getDCELPointCount(dcel);
            setupGTKPolyPoints(dcelVerticesCount,
                (void *) dcel,
                (double (*)(void *, int)) &getDCELVertexX,
                (double (*)(void *, int)) &getDCELVertexY
            );

            /* Visualisation: DCEL lines */
            int dcelEdgeCount = getDCELEdgeCount(dcel);
            setupGTKPolyLines(dcelEdgeCount,
                (void *) dcel,
                (int (*)(void *, int)) &getDCELEdgeVertexStart,
                (int (*)(void *, int)) &getDCELEdgeVertexPairStart,
                (int (*)(void *, int)) &getDCELEdgeVertexEnd,
                (int (*)(void *, int)) &getDCELEdgeVertexPairEnd,
                (int (*)(void *, int)) &DCELhasEdge,
                (int (*)(void *, int)) &DCELhasEdgePair
            );
            #endif
            
            break;
            
        case STAGE_1:
            pointpairfile = argv[2];
            outputfile = argv[3];
            
            ppFile = fopen(pointpairfile, "r");
            assert(ppFile);
            outFile = fopen(outputfile, "w");
            assert(outFile);
            
            while((currentBisector = readNextBisector(ppFile)) != NULL){
                char *bisectorString = getBisectorEquation(currentBisector);
                if (bisectorString){
                    fprintf(outFile, "%s\n", bisectorString);
                }
                
                freeBisector(currentBisector);
                if (bisectorString){
                    free(bisectorString);
                }
            }
            
            fclose(ppFile);
            fclose(outFile);
            break;
            
        case STAGE_2:
            pointpairfile = argv[2];
            polygonfile = argv[3];
            outputfile = argv[4];
            
            ppFile = fopen(pointpairfile, "r");
            assert(ppFile);
            outFile = fopen(outputfile, "w");
            assert(outFile);
            
            dcel = readPolygonFile(polygonfile);
            
            while((currentBisector = readNextBisector(ppFile)) != NULL){
                struct intersection *intersection = 
                    getIntersection(currentBisector, dcel, DEFAULT_FACE, DEFAULTMINLENGTH);
                char *intersectionString = getIntersectionString(intersection);
                fprintf(outFile, "%s\n", intersectionString);
                
                freeIntersection(intersection);
                free(intersectionString);
                freeBisector(currentBisector);
            }
            
            fclose(ppFile);
            fclose(outFile);
            
            break;
            
        case STAGE_3:
            datafile = argv[2];
            polygonfile = argv[3];
            outputfile = argv[4];
            
            watchTowers = readDataFile(datafile, &watchTowerCount);
            dcel = readPolygonFile(polygonfile);
            for(int i = 0; i < watchTowerCount; i++){
                // We add one watchTower at a time.
                incrementalVoronoi(dcel, watchTowers[i]);
            }
            
            #ifdef USE_GUI
            /* Visualisation: Points */
            pointCount = watchTowerCount;
            setupGTKPoints(pointCount, 
                (void *) watchTowers,
                (char *(*)(void *, int)) &getWTDataString, 
                (void (*)(char *)) &freeWTDataString,
                (double (*)(int, void *)) &getWTx,
                (double (*)(int, void *)) &getWTy
            );

            /* Visualisation: DCEL points */
            dcelVerticesCount = getDCELPointCount(dcel);
            setupGTKPolyPoints(dcelVerticesCount,
                (void *) dcel,
                (double (*)(void *, int)) &getDCELVertexX,
                (double (*)(void *, int)) &getDCELVertexY
            );

            /* Visualisation: DCEL lines */
            dcelEdgeCount = getDCELEdgeCount(dcel);
            setupGTKPolyLines(dcelEdgeCount,
                (void *) dcel,
                (int (*)(void *, int)) &getDCELEdgeVertexStart,
                (int (*)(void *, int)) &getDCELEdgeVertexPairStart,
                (int (*)(void *, int)) &getDCELEdgeVertexEnd,
                (int (*)(void *, int)) &getDCELEdgeVertexPairEnd,
                (int (*)(void *, int)) &DCELhasEdge,
                (int (*)(void *, int)) &DCELhasEdgePair
            );
            #endif
            
            break;
            
        case STAGE_4:
            datafile = argv[2];
            polygonfile = argv[3];
            outputfile = argv[4];
            
            watchTowers = readDataFile(datafile, &watchTowerCount);
            dcel = readPolygonFile(polygonfile);
            for(int i = 0; i < watchTowerCount; i++){
                // We add one watchTower at a time.

                incrementalVoronoi(dcel, watchTowers[i]);
            }
            
            #ifdef USE_GUI
            /* Visualisation: Points */
            pointCount = watchTowerCount;
            setupGTKPoints(pointCount, 
                (void *) watchTowers,
                (char *(*)(void *, int)) &getWTDataString, 
                (void (*)(char *)) &freeWTDataString,
                (double (*)(int, void *)) &getWTx,
                (double (*)(int, void *)) &getWTy
            );

            /* Visualisation: DCEL points */
            dcelVerticesCount = getDCELPointCount(dcel);
            setupGTKPolyPoints(dcelVerticesCount,
                (void *) dcel,
                (double (*)(void *, int)) &getDCELVertexX,
                (double (*)(void *, int)) &getDCELVertexY
            );

            /* Visualisation: DCEL lines */
            dcelEdgeCount = getDCELEdgeCount(dcel);
            setupGTKPolyLines(dcelEdgeCount,
                (void *) dcel,
                (int (*)(void *, int)) &getDCELEdgeVertexStart,
                (int (*)(void *, int)) &getDCELEdgeVertexPairStart,
                (int (*)(void *, int)) &getDCELEdgeVertexEnd,
                (int (*)(void *, int)) &getDCELEdgeVertexPairEnd,
                (int (*)(void *, int)) &DCELhasEdge,
                (int (*)(void *, int)) &DCELhasEdgePair
            );
            #endif
            
            break;
            
        case STAGE_ERROR:
            fprintf(stderr, "Failed to read arguments");
            exit(EXIT_FAILURE);
            break;
        default:
            // UNSET & STAGE_ERROR
            fprintf(stderr, "Unhandled stage (%s).\n", getStageString(stage));
    }
    
    
    
    #ifdef USE_GUI
    /* Set click action and data. */
    /* Show location of mouse-click */
    /* setupGTKAction((void (*)(void *, double, double)) &doAction, NULL); */
    int gtkActionData[3];
    // gtkActionData[0] = WT_CHANGE;
    // gtkActionData[0] = POLY_CHANGE;
    gtkActionData[0] = EDGE_CHANGE;
    gtkActionData[1] = 0;
    
    setupGTKAction((void (*)(void *, double, double)) &doCyleOnClick, 
                   (void *) &(gtkActionData[0]));
    
    runGTKInteraction();
    /* NB: GTK has memory leaks so no need to have clean output. */
    freeGTK();
    #endif
    
    switch(stage){
        case STAGE_0:
            outputResult(outputfile, watchTowers, watchTowerCount, dcel);
            freeWatchTowers(watchTowers, watchTowerCount);
            freeDCEL(dcel);
            break;
            
        case STAGE_1:
            break;
            
        case STAGE_2:
            freeDCEL(dcel);
            break;
            
        case STAGE_3:
            outputResultDiameter(outputfile, watchTowers, watchTowerCount, dcel);
            freeWatchTowers(watchTowers, watchTowerCount);
            freeDCEL(dcel);
            break;
            
        case STAGE_4:
            outputResultDiameterSorted(outputfile, watchTowers, watchTowerCount, dcel);
            freeWatchTowers(watchTowers, watchTowerCount);
            freeDCEL(dcel);
            break;
            
        default:
            // UNSET & STAGE_ERROR
            fprintf(stderr, "Unexpected unhandled stage at output/free step (%s).\n", getStageString(stage));       
    }
    
    return EXIT_SUCCESS;
}

