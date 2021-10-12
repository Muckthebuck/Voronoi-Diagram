#include "watchtowerStruct.h"

/* Reads the given csv datafile, allocating memory for an array of pointers, storing 
the number of read values into watchTowerCount. */
struct watchtowerStruct **readDataFile(char *datafileName, int *watchTowerCount);

/* Frees all the memory allocated by the readDataFile function. */
void freeWatchTowers(struct watchtowerStruct **wts, int wtCount);
