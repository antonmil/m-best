#ifndef PARTICLES_NEW2_H
#define PARTICLES_NEW2_H

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "rngs.h"
#include <string.h>

typedef struct {
    int nMembers;        // number of cluster members;
    int *trajectory;
    int *tempProb;
    double *trajSum;        // transition probability
    double *cumTranProb;    // cumulative transition probability
    double sumAffinity;
    bool termFlag;
    int length;
} cluster;

typedef struct {
    int nMembers;
    int *members;
} groupList;

typedef struct {
	groupList *groupList1;
	groupList *groupList2;
	int *groupRef;			// nBins by 2. what group the match belongs
	int nGroup1;
	int nGroup2;
} confGroupInfo;

typedef struct LLNode{
	int data;
	LLNode *next;
} LLNode;

void initDistribution(int **xPfInit);

void group(int *particleId, cluster** newClusters);
void group2(int *particleId, cluster** newClusters, int *prevClusterId, int *nextParticle, int *tempTraj, int *sortParticleId);
void setCluster(int next, int nMembers, cluster *newCluster); //, cluster *prevCluster);
void setCluster2(int next, int nMembers, cluster *newCluster, cluster *prevCluster, int *tempTraj);
void setGroupList(int *tempGroup1, int *tempGroup2);
void freeClusters(cluster* clusters, int nClusters);
void freeGroupList();

void desample(cluster* newClusters, cluster* newClustersDD);
void desetCluster(int theMember, int theMemberId, double minSum, int i, cluster *newClusters, cluster *newClustersDD);
//int makeClusterId(int* clusterId, cluster* newClusters, cluster** saveCluster, double* nTrajectory, double* trajectoryIdx, int* trajCount, int* nTrajCount);
int makeClusterId(int* clusterId, cluster* newClusters, cluster** saveCluster, double* nTrajectory, double* trajectoryIdx, int* trajCount, int* nTrajCount);

int sample(double* cumProb, int n);
void normalizeWeights(double *x, int nX);
int countUnique(int *x, int nX, int **nMembers);
int compare (const void * a, const void * b);
bool isEqualComb(int *comb1, int *comb2);
int* radixSort(int* tempTraj);
void insertNew(int xPfPred, int tempClusterId, int* tempTraj, int j);
double kernel(double x);
double getTemperT();

void LLInsert(LLNode **bucketEnds, int bucketIdx, int data);
void freeLL(LLNode *bucket);

#endif