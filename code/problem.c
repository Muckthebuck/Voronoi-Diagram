#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "problem.h"

// ./voronoi2 0 datafile polygonfile outputfile < splitsfile
#define STAGE_0_ARG_COUNT 5
// ./voronoi2 1 pointpairs outputfile
#define STAGE_1_ARG_COUNT 4
// ./voronoi2 2 pointpairs polygonfile outputfile
#define STAGE_2_ARG_COUNT 5
// ./voronoi2 3 datafile polygonfile outputfile
#define STAGE_3_ARG_COUNT 5
// ./voronoi2 4 datafile polygonfile outputfile
#define STAGE_4_ARG_COUNT 5
// Shouldn't be needed.
#define UNSET_ARG_COUNT 2

enum assignmentStage getArgumentMode(char *mode){
    /*
        STAGE_0 = 0,
        STAGE_1 = 1,
        STAGE_2 = 2,
        STAGE_3 = 3,
        STAGE_4 = 4,
        STAGE_ERROR = (-1)
    */
    
    if(strcmp(mode, "0") == 0){
        return STAGE_0;
    } else if (strcmp(mode, "1") == 0){
        return STAGE_1;
    } else if (strcmp(mode, "2") == 0){
        return STAGE_2;
    } else if (strcmp(mode, "3") == 0){
        return STAGE_3;
    } else if (strcmp(mode, "4") == 0){
        return STAGE_4;
    } else {
        return STAGE_ERROR;
    }
}

char *getStageString(enum assignmentStage stage){
    switch(stage){
        case STAGE_0:
            return "STAGE_0";
            break;
        case STAGE_1:
            return "STAGE_1";
            break;
        case STAGE_2:
            return "STAGE_2";
            break;
        case STAGE_3:
            return "STAGE_3";
            break;
        case STAGE_4:
            return "STAGE_4";
            break;
        case STAGE_ERROR:
            return "STAGE_ERROR";
            break;
        case UNSET:
            return "UNSET";
            break;
        default:
            fprintf(stderr, "Stage was not handled ENUM value (%d)\n", stage);
            exit(EXIT_FAILURE);
            return "ERROR";
    }
}

void printArgError(char *error, int argCount, enum assignmentStage stage){
        fprintf(stderr, "%s, run with:\n\n", error);
        
        if(stage == UNSET || stage == STAGE_ERROR || stage == STAGE_0){
            /* Mode 0: Voronoi1 */
            fprintf(stderr, "\t./voronoi2 0 datafile polygonfile outputfile < splitsfile\n"
                            "\t\tTo generate the polygon in polygonfile, layout the points in \n"
                            "\t\tthe datafile on the polygon, apply the splits in the splitsfile \n"
                            "\t\tand output to the outputfile\n\n");
        }
        if(stage == UNSET || stage == STAGE_ERROR || stage == STAGE_1){
            /* Mode 1: Voronoi2 Stage 1 - Print bisector equations */
            fprintf(stderr, "\t./voronoi2 1 pointpairs outputfile\n"
                            "\t\tTo output the bisectors for the given pointpairs to the given\n"
                            "\t\toutputfile.\n\n");
        }
        if(stage == UNSET || stage == STAGE_ERROR || stage == STAGE_2){
            /* Mode 2: Voronoi2 Stage 2 - Print intersections between bisectors and polygon. */
            fprintf(stderr, "\t./voronoi2 2 pointpairs polygonfile outputfile\n"
                            "\t\tTo output the intersections the given bisectors make with the polygon \n"
                            "\t\tto the given outputfile.\n\n");
        }
        if(stage == UNSET || stage == STAGE_ERROR || stage == STAGE_3){
            /* Mode 3: Voronoi2 Stage 3 - Print diameter of each face, constructing using incremental method. */
            fprintf(stderr, "\t./voronoi2 3 datafile polygonfile outputfile\n"
                            "\t\tTo generate the initial polygon in polygonfile, layout the watchtower \n"
                            "\t\tpoints in the datafile on the polygon, generate the voronoi diagram, adding\n"
                            "\t\tmodifying the polygon so that a voronoi diagram is constructed, seperating \n"
                            "\t\tall points into their own cell with all points in the cell being closest to\n"
                            "\t\tthe watchtower in the cell. The watchtower data and its diameter is\n"
                            "\t\toutput to the outputfile\n\n");
        }
        if(stage == UNSET || stage == STAGE_ERROR || stage == STAGE_4){
            /* Mode 4: Voronoi2 Stage 4 - Sort by diameter. */
            fprintf(stderr, "\t./voronoi2 4 datafile polygonfile outputfile\n"
                            "\t\tSame as stage 3, but output watchtowers are in order of diameter, smallest \n"
                            "\t\tto largest.\n\n");
        }
        exit(EXIT_FAILURE);
}

int expectedArgs(enum assignmentStage stage){
    switch(stage){
        case STAGE_0:
            return STAGE_0_ARG_COUNT;
        case STAGE_1:
            return STAGE_1_ARG_COUNT;
        case STAGE_2:
            return STAGE_2_ARG_COUNT;
        case STAGE_3:
            return STAGE_3_ARG_COUNT;
        case STAGE_4:
            return STAGE_4_ARG_COUNT;
        case STAGE_ERROR:
            return (-1);
        case UNSET:
            return UNSET_ARG_COUNT;
        default:
            fprintf(stderr, "Stage was not handled ENUM value (%d)\n", stage);
            exit(EXIT_FAILURE);
            return (-1);
    }
}

int checkArgCount(enum assignmentStage stage, int argCount){
    if(argCount >= expectedArgs(stage)){
        return 1;
    } else {
        return 0;
    }
}
