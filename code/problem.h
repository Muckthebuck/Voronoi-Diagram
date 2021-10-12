enum assignmentStage;


#ifndef PROBLEMSTAGE_ENUM
#define PROBLEMSTAGE_ENUM
/* We want this to be visible to anything using this compilation object. */
enum assignmentStage {
    STAGE_0 = 0,
    STAGE_1 = 1,
    STAGE_2 = 2,
    STAGE_3 = 3,
    STAGE_4 = 4,
    STAGE_ERROR = (-1),
    UNSET = (-2)
};
#endif

/* 
    Reads the string and interprets the string as the correct mode if
    possible, returning this stage. Returns STAGE_ERROR if the argument
    is not a valid stage string.
*/
enum assignmentStage getArgumentMode(char *mode);

/* 
    Gets a string representing the stage.
*/
char *getStageString(enum assignmentStage stage);

/*
    Prints an error and the usage expected.
*/
void printArgError(char *error, int argCount, enum assignmentStage stage);

/*
    Checks if the argument count is correct for the given stage.
    Returns 1 if the argument count is correct, 0 otherwise.
*/
int checkArgCount(enum assignmentStage stage, int argCount);
