#include <stdio.h>
#include <stdlib.h>
#include "RCA.h"
#include <unistd.h>
#include <sys/param.h>
#include <string.h>
#include <math.h>

int * mMembership;
int mMergeDim[2];
int ** mMerges; 
double ** mDataset;
char * mFilename;

void RCA_R(double ** data_set, int numObs, int numVars, int **membership, double * mod, double * returnStats, int *merge_dim, int *** merges, int bootstrap, double z_score) {
    *membership=NULL;
    double modularity;
    double stats[2];
    static const char* c_assignments = "member";
    static const char* c_merges = "merge";
    char buf[50];

    RCA(data_set, numObs, numVars, membership, &modularity, merges, merge_dim, stats, bootstrap, z_score);

    int group_num = 0;

    for (int i=0; i<numObs; i++) {
        if(group_num<(*membership)[i]) {
            group_num = (*membership)[i];
        }
    }
    group_num++;

    *mod=modularity;
    returnStats[0]=stats[0];
    returnStats[1]=stats[1];
}

void RCA_R_wrapper(double * inputData1D, int * numObs, int * numVars, int * bootstrap, double * z_score, int *mem, int * mergeDim, double * mod, double * stats){
    /*Initialize member matrix.*/
	mMerges=(int**)malloc(MAX_COL*sizeof(int*));
    for(int i=0;i<MAX_COL;i++){
        *mMerges=(int*)malloc(MAX_COL*sizeof(int));
    }
	/*Change input data from 1D to 2D.*/
	double ** dataset2D=(double**)malloc((*numObs)*sizeof(double*));
	for(int i=0;i<(*numObs);i++){
		dataset2D[i]=(double*)malloc((*numVars)*sizeof(double));
	}
	for(int i=0;i<(*numObs);i++){
		for(int j=0;j<(*numVars);j++){
			dataset2D[i][j]=inputData1D[i*(*numVars)+j];
		}
	}
    RCA_R(dataset2D, (*numObs), (*numVars), &mem, mod, stats, mergeDim, &mMerges, *bootstrap, *z_score);
    mMembership=mem;
	mMergeDim[0]=mergeDim[0];
    mMergeDim[1]=mergeDim[1];
}

/* Get values using a getter function because members array is too big. 
*/
void get(int * ind, int * val){
    *val = mMembership[*ind];
}

/* Get values using a getter function because: dimension of merges is unkown. */
void getMerges(int * merges1D){
    int row = mMergeDim[0];
    int col = mMergeDim[1];
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
            merges1D[j*row+i]=mMerges[i][j];
        }
    }
}
