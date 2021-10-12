#include "output.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "watchtowerStruct.c"
#include "dcel.h"

#define NOFACE (-1)

typedef struct {
	double distance;
	char* watchtower_id;
	int idx;
}task4_t;

void outputResult(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){
    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    int i, j;
    int population = 0;
    if(!dcel){
        /* Simple path - avoids DCEL entirely. */
        for(i = 0; i < wtCount; i++){
            if(! wts[i]){
                continue;
            }
            fprintf(outputfile,
                "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
                "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
                (wts[i])->watchtowerID, 
                (wts[i])->postcode, 
                (wts[i])->populationServed, 
                (wts[i])->contact, 
                (wts[i])->x, 
                (wts[i])->y);
            population += wts[i]->populationServed;
        }
        fprintf(outputfile,
            "Face undefined population served: %d\n", population);
    } else {
        int faceCount = getFaceCount(dcel);
        int *faceMembership = (int *) malloc(sizeof(int)*wtCount);
        assert(faceMembership);
        for(i = 0; i < wtCount; i++){
            faceMembership[i] = NOFACE;
            for(j = 0; j < faceCount; j++){
                if(inFace(dcel, wts[i]->x, wts[i]->y, j)){
                    // if(faceMembership[i] != NOFACE){
                    //     printf("WT %d on edge of %d and %d!\n", i, faceMembership[i], j);
                    // }
                    faceMembership[i] = j;
                    break;
                }
            }
        }
        int *faceSize = (int *) malloc(sizeof(int)*(faceCount + 1));
        assert(faceSize);
        for(i = 0; i < (faceCount + 1); i++){
            faceSize[i] = 0;
        }
        for(i = 0; i < wtCount; i++){
            if(faceMembership[i] == NOFACE){
                faceSize[faceCount] = faceSize[faceCount] + 1;
            } else {
                faceSize[faceMembership[i]] = faceSize[faceMembership[i]] + 1;
            }
        }
        int **faceData = (int **) malloc(sizeof(int *)*(faceCount + 1));
        assert(faceData);
        for(i = 0; i < (faceCount + 1); i++){
            faceData[i] = (int *) malloc(sizeof(int)*faceSize[i]);
            assert(faceData[i] || faceSize[i] == 0);
            /* Reset the size of arrays so we can put them back in. */
            faceSize[i] = 0;
        }
        /* Put them into arrays. */
        for(i = 0; i < wtCount; i++){
            if(faceMembership[i] == NOFACE){
                faceData[faceCount][faceSize[faceCount]] = i;
                faceSize[faceCount] = faceSize[faceCount] + 1;
            } else {
                faceData[faceMembership[i]][faceSize[faceMembership[i]]] = i;
                faceSize[faceMembership[i]] = faceSize[faceMembership[i]] + 1;
            }
        }
        if(faceMembership){
            free(faceMembership);
        }
        /* Output */
        int *populationTotals = (int *) malloc(sizeof(int) * (faceCount + 1));
        assert(populationTotals);
        for(i = 0; i < (faceCount + 1); i++){
            populationTotals[i] = 0;
            // if(i != faceCount || faceSize[i] != 0){
            if(i != faceCount){
                fprintf(outputfile, "%d\n", i);
                for(j = 0; j < faceSize[i]; j++){
                    fprintf(outputfile,
                        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
                        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
                        (wts[faceData[i][j]])->watchtowerID, 
                        (wts[faceData[i][j]])->postcode, 
                        (wts[faceData[i][j]])->populationServed, 
                        (wts[faceData[i][j]])->contact, 
                        (wts[faceData[i][j]])->x, 
                        (wts[faceData[i][j]])->y);
                    populationTotals[i] += wts[faceData[i][j]]->populationServed;
                }
            }
        }
        
        for(i = 0; i < (faceCount); i++){
            fprintf(outputfile,
                "Face %d population served: %d\n", i, populationTotals[i]);
        }
        // if(faceSize[faceCount] > 0){
        //     fprintf(outputfile,
        //         "Face undefined population served: %d\n", populationTotals[faceCount]);
        // }
        if(populationTotals){
            free(populationTotals);
        }
        if(faceData){
            for(i = 0; i < (faceCount + 1); i++){
                if(faceData[i]){
                    free(faceData[i]);
                }
            }
            free(faceData);
        }
        if(faceSize){
            free(faceSize);
        }
    }
    fclose(outputfile);
}

void outputResultDiameter(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){
    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    int i;
    
    /* Must have DCEL. */
    assert(dcel);
    
    for(i = 0; i < wtCount; i++){
        if(! wts[i]){
            continue;
        }
        double diameter = getDiameter(dcel, wts[i]->face);
        fprintf(outputfile,
            "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
            "Watchtower Point of Contact Name: %s, x: %f, y: %f, Diameter of Cell: %f\n",
            (wts[i])->watchtowerID, 
            (wts[i])->postcode, 
            (wts[i])->populationServed, 
            (wts[i])->contact, 
            (wts[i])->x, 
            (wts[i])->y,
            diameter);
    }
    
    fclose(outputfile);
}

void outputResultDiameterSorted(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){
    /* TODO: 
        Fill in for Task 4:
        Make this sorted using insertion sort.
    */
    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    int i;
    
    /* Must have DCEL. */
    assert(dcel);
    task4_t arr[wtCount];
    for(i = 0; i < wtCount; i++){
    	if(! wts[i]){
    		continue;
    	}
    	double diameter = getDiameter(dcel, wts[i]->face);
    	arr[i].distance = diameter;
    	arr[i].watchtower_id = wts[i]->watchtowerID;
    	arr[i].idx = i;
    	int k, j;
    	for (k =1;k<i+1;k++){
			double v1 = arr[k].distance;
			char*  v2 = arr[k].watchtower_id;
			int v3 = arr[k].idx;
			j = k-1;
			while(j>=0 && ((arr[j].distance>v1) ||
			     ((arr[j].distance==v1)&&arr[j].watchtower_id>v2))){
				arr[j+1].distance = arr[j].distance;
				arr[j+1].watchtower_id = arr[j].watchtower_id;
				arr[j+1].idx = arr[j].idx;
				j=j-1;
			}
			arr[j+1].distance=v1;
			arr[j+1].watchtower_id= v2;
			arr[j+1].idx = v3;
    	}
    }
    for(i = 0; i < wtCount; i++){
        if(! wts[i]){
            continue;
        }

        fprintf(outputfile,
            "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
            "Watchtower Point of Contact Name: %s, x: %f, y: %f, Diameter of Cell: %f\n",
            (wts[arr[i].idx])->watchtowerID,
            (wts[arr[i].idx])->postcode,
            (wts[arr[i].idx])->populationServed,
            (wts[arr[i].idx])->contact,
            (wts[arr[i].idx])->x,
            (wts[arr[i].idx])->y,
            arr[i].distance);
    }
    
    fclose(outputfile);
}

char *getWTDataString(struct watchtowerStruct **wts, int wtIndex){
    char *s = NULL;
    /* Get size */
    size_t requiredSpace = snprintf(NULL, 0, 
        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
        (wts[wtIndex])->watchtowerID, 
        (wts[wtIndex])->postcode, 
        (wts[wtIndex])->populationServed, 
        (wts[wtIndex])->contact, 
        (wts[wtIndex])->x, 
        (wts[wtIndex])->y);
    s = (char *) malloc(sizeof(char) * (requiredSpace + 1));
    assert(s);
    sprintf(s, 
        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
        (wts[wtIndex])->watchtowerID, 
        (wts[wtIndex])->postcode, 
        (wts[wtIndex])->populationServed, 
        (wts[wtIndex])->contact, 
        (wts[wtIndex])->x, 
        (wts[wtIndex])->y);
    
    return s;
}

void freeWTDataString(char *s){
    if(s){
        free(s);
    }
}

double getWTx(int index, struct watchtowerStruct **wts){
    if(!wts || !(wts[index])){
        return 0;
    }
    
    return (wts[index])->x;
}

double getWTy(int index, struct watchtowerStruct **wts){
    if(!wts || !(wts[index])){
        return 0;
    }
    
    return (wts[index])->y;
}
