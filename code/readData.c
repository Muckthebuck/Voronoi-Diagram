#include "watchtowerStruct.h"
// Needs internals of struct to work
#include "watchtowerStruct.c"
#include "readData.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define INITIALWTS 1

#define FINAL_FIELD 5

struct watchtowerStruct **readDataFile(char *datafileName, int *watchTowerCount){
    assert(datafileName);
    FILE *file = fopen(datafileName, "r");
    assert(file);
    char *line = NULL;
    int progress;
    size_t lineAlloced = 0;
    int first = 1;
    int fieldStart;
    int fieldProgress;
    
    struct watchtowerStruct **wts = NULL;
    int wtsAlloced = 0;
    int wtsUsed = 0;
    
    struct watchtowerStruct *currentWT = NULL;
    
    while(getline(&line, &lineAlloced, file) != (-1)){
        /* Ignore header for simplicity. */
        if(first){
            first = 0;
            continue;
        }
        
        /* Remove new line character if present */
        if(line[strlen(line)] == '\n'){
            line[strlen(line)] = '\0';
        }
        
        fieldStart = 0;
        progress = 0;
        fieldProgress = 0;
        
        /* Ensure we have space for the new watchtower. */
        if(wtsAlloced == 0){
            wts = (struct watchtowerStruct **) 
                malloc(sizeof(struct watchtowerStruct *) * INITIALWTS);
            assert(wts);
            wtsAlloced = INITIALWTS;
        } else if ((wtsUsed + 1) > wtsAlloced){
            wts = (struct watchtowerStruct **) 
                realloc(wts, sizeof(struct watchtowerStruct *)*wtsAlloced*2);
            assert(wts);
            wtsAlloced *= 2;
        }
        currentWT = (struct watchtowerStruct *) malloc(sizeof(struct watchtowerStruct));
        assert(currentWT);
        
        while(line[progress] != '\0'){
            
            if(line[progress] == ','){
                // Field ended, extract field
                switch(fieldProgress){
                    case 0:
                        line[progress] = '\0';
                        //currentWT->watchTowerID = (char *) 
                        //    malloc((strlen(line + fieldStart) + 1) * sizeof(char));
                        // assert(currentWT->watchTowerID);
                        // strncpy(currentWT->watchTowerID, line + fieldStart, 
                        //    strlen(line + fieldStart) + 1);
                        currentWT->watchtowerID = strdup(line + fieldStart);
                        
                        fieldStart = progress + 1;
                        fieldProgress++;
                        break;
                        
                    case 1:
                        line[progress] = '\0';
                        currentWT->postcode = strdup(line + fieldStart);
                        fieldStart = progress + 1;
                        fieldProgress++;
                        break;
                        
                    case 2:
                        line[progress] = '\0';
                        assert(sscanf(line + fieldStart, "%d", 
                            &(currentWT->populationServed)) == 1);
                        fieldStart = progress + 1;
                        fieldProgress++;
                    break;
                    
                    case 3:
                        line[progress] = '\0';
                        currentWT->contact = strdup(line + fieldStart);
                        fieldStart = progress + 1;
                        fieldProgress++;
                        break;
                    break;
                    
                    case 4:
                        line[progress] = '\0';
                        assert(sscanf(line + fieldStart, "%lf", 
                            &(currentWT->x)) == 1);
                        fieldStart = progress + 1;
                        fieldProgress++;
                    break;
                    
                    //case 5:
                    // Handled by '\0' case.
                    //break;
                        
                    default:
                        fprintf(stderr, "Error: too many fields in line %s\n",
                            currentWT->watchtowerID);
                }
            }
            progress++;
        }
        if(line[progress] == '\0'){
            assert(fieldProgress == FINAL_FIELD);
            assert(sscanf(line + fieldStart, "%lf", 
                &(currentWT->y)) == 1);
            fieldStart = progress + 1;
            fieldProgress++;
        }
        
        wts[wtsUsed] = currentWT;
        wtsUsed++;
    }
    
    if(line){
        free(line);
    }
    fclose(file);
    
    *watchTowerCount = wtsUsed;
    return wts;
}

void freeWatchTowers(struct watchtowerStruct **wts, int wtCount){
    if(! wts){
        return;
    }
    int i;
    for(i = 0; i < wtCount; i++){
        if(wts[i]){
            if(wts[i]->watchtowerID){
                free(wts[i]->watchtowerID);
            }
            if(wts[i]->postcode){
                free(wts[i]->postcode);
            }
            if(wts[i]->contact){
                free(wts[i]->contact);
            }
            free(wts[i]);
            wts[i] = NULL;
        }
    }
    free(wts);
}
