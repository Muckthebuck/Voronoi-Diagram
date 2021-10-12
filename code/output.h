#include "watchtowerStruct.h"
#include "dcel.h"

/* Outputs all data to file with given filename in specification format. */
void outputResult(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel);

/* Outputs all data to file with given filename in specification format with diameter. */
void outputResultDiameter(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel);

/* Outputs all data to file with given filename in specification format, sorted by diameter. */
void outputResultDiameterSorted(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel);

/* Returns a string for the given watchtower. */
char *getWTDataString(struct watchtowerStruct **wts, int wtIndex);

/* Simple free string wrapper. */
void freeWTDataString(char *s);

/* Get the x value from the wts array given. */
double getWTx(int index, struct watchtowerStruct **wts);

/* Get the y value from the wts array given. */
double getWTy(int index, struct watchtowerStruct **wts);
