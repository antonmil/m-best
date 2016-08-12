#include "particles_eccv.h"
///////////////////////////////////////////////////////////////////////////
                        /* Function Definitions */
///////////////////////////////////////////////////////////////////////////

extern confGroupInfo confInfo;
extern int nCandMatch;
extern double *M;
extern int nParticle;
extern int l;
extern cluster *clusters;
extern int nCluster;
extern int nNewCluster;
extern int *rankM;
extern int nWindow;
extern double temperature;

/* initialize distribution */
void initDistribution(int **xPfInit) {
    *xPfInit = (int*)malloc(nParticle * sizeof(int));
	double *expM = (double*)malloc(nCandMatch*nCandMatch * sizeof(double));
    double *initProb = (double*)malloc(nCandMatch * sizeof(double));
	memset(initProb, 0, nCandMatch * sizeof(double));
    double *cumInitProb = (double*)malloc(nCandMatch * sizeof(double));
    double sumM = 0;
    
    int i, j, i1, i2, j1, j2;

    for(i = 0; i < nCandMatch; i++) {
		int nCandMatchI = nCandMatch*i;
        for(j = 0; j < nCandMatch; j++) {
            expM[nCandMatchI + j] = exp(M[nCandMatchI + j]/getTemperT());
        }
    }
    
    for(i = 0; i < confInfo.nGroup1; i++) {
        for(i1 = 0; i1 < confInfo.groupList1[i].nMembers; i1++) {
            for(i2 = 0; i2 < confInfo.groupList1[i].nMembers; i2++) {
                expM[ confInfo.groupList1[i].members[i1] + nCandMatch*confInfo.groupList1[i].members[i2] ] = 0;
            }
        }
    }
    for(j = 0; j < confInfo.nGroup2; j++) {
        for(j1 = 0; j1 < confInfo.groupList2[j].nMembers; j1++) {
            for(j2 = 0; j2 < confInfo.groupList2[j].nMembers; j2++) {
                expM[ confInfo.groupList2[j].members[j1] + nCandMatch*confInfo.groupList2[j].members[j2] ] = 0;
            }
        }
    }
    
    // get initProb from M
    for(i = 0; i < nCandMatch; i++) {
		int nCandMatchI = nCandMatch*i;
        for(j = 0; j < nCandMatch; j++) {
//             initProb[i] += M[i + nCandMatch*j];
//             initProb[i] += exp(M[i + nCandMatch*j]) * M[i + nCandMatch*j];
            initProb[i] += expM[nCandMatchI + j];
        }
        sumM += initProb[i];
    }
    free(expM);
    
    // normalize
    for(i = 0; i < nCandMatch; i++) {
        initProb[i] = initProb[i] / sumM;
    }
    
    // get cumInitProb from initProb
    cumInitProb[0] = initProb[0];
    for(i = 1; i < nCandMatch; i++) {
        cumInitProb[i] = cumInitProb[i-1] + initProb[i];
    }
        
    // sample
    for(i = 0; i < nParticle ; i++) {
        *(*xPfInit+i) = sample(cumInitProb, nCandMatch);
    }
    free(initProb);
    free(cumInitProb);
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

/* group clusters */
void group(int *particleId, cluster** newClusters) {

	int i, j;
    int *nMembers;
	int tempNext;
    
    // get unique cluster information
	qsort(particleId, nParticle, sizeof(int), compare);
	nNewCluster = countUnique(particleId, nParticle, &nMembers);
        
    // make new clusters
    *newClusters = (cluster*)malloc(nNewCluster * sizeof(cluster));
    tempNext = particleId[0] % nCandMatch;
    
    j = 0;
    setCluster(tempNext, nMembers[j], *newClusters) ;
    for(i = 1; i < nParticle; i++) {
        if(particleId[i] != particleId[i-1]) {
            j++;
            tempNext = particleId[i] % nCandMatch;
            setCluster(tempNext, nMembers[j], (*newClusters)+j );
        }
    }
}

/* group clusters */
void group2(int *outIndex, cluster** newClusters, int *prevClusterId, int *nextParticle, int *tempTraj, int *sortParticleId) {

	int i, j, k;
	int *nMembers;
	int *tempComb;
    
    // get unique cluster information
    qsort(outIndex, nParticle, sizeof(int), compare);
    nNewCluster = countUnique(outIndex, nParticle, &nMembers);
        
    // make new clusters
    *newClusters = (cluster*)malloc(nNewCluster * sizeof(cluster));
    
    j = 0;
    tempComb = (int*)malloc(l * sizeof(int));
    for(k = 0 ; k < l; k++) {
        tempComb[k] = tempTraj[sortParticleId[outIndex[j]] + k*nParticle];
    }

    setCluster2(nextParticle[outIndex[0]], nMembers[0], *newClusters, (clusters + prevClusterId[ outIndex[0] ]), tempComb);
    free(tempComb);
    for(i = 1; i < nParticle; i++) {
        if(outIndex[i] != outIndex[i-1]) {
            j++;
            tempComb = (int*)malloc(l * sizeof(int));
            for(k = 0 ; k < l; k++) {
                tempComb[k] = tempTraj[sortParticleId[outIndex[i]] + k*nParticle];
            }

            setCluster2(nextParticle[outIndex[i]], nMembers[j], (*newClusters)+j, (clusters + prevClusterId[ outIndex[i] ]), tempComb);
            free(tempComb);
        }
    }
}

/* set cluster */
void setCluster(int next, int nMembers, cluster *newCluster ) { //, cluster *prevCluster) {

	double *MNext = (double*)malloc(nCandMatch *sizeof(double));
    double *tranProb;
    int tempProbSum = 0;
	double sum = 0;
	int tempLastId;
    
	int i, j;
	
    for(i = 0; i < nCandMatch; i++) {
        MNext[i] = M[next + nCandMatch*i];
    }
        
    // set the number of members
    (*newCluster).nMembers = nMembers;
    
    // set the trajectory
    (*newCluster).trajectory = (int*)malloc(sizeof(int));
    (*newCluster).trajectory[0] = next;
    
    // set the sumAffinity
    (*newCluster).sumAffinity = 0;

    // set the transition probability
    (*newCluster).trajSum = (double*)malloc(nCandMatch * sizeof(double));
    for(i = 0; i < nCandMatch; i++) {
        (*newCluster).trajSum[i] = MNext[i];
    }
    
    // set the tempProb
    (*newCluster).tempProb = (int*)malloc(nCandMatch * sizeof(int));
    for(i = 0; i < nCandMatch; i++) {
        (*newCluster).tempProb[i] = 1;
    }
    for(i = 0; i < confInfo.groupList1[confInfo.groupRef[next]].nMembers; i++) {
        (*newCluster).tempProb[ confInfo.groupList1[confInfo.groupRef[next]].members[i] ] = 0;
    }
    for(i = 0; i < confInfo.groupList2[confInfo.groupRef[next + nCandMatch]].nMembers; i++) {
        (*newCluster).tempProb[ confInfo.groupList2[confInfo.groupRef[next + nCandMatch]].members[i] ] = 0;
    }
    
    (*newCluster).termFlag = false;
    for(i = 0;i < nCandMatch; i++) {
        tempProbSum += (*newCluster).tempProb[i];
    }
    if(tempProbSum == 0)    (*newCluster).termFlag = true;
  
    tranProb = (double*)malloc(nCandMatch * sizeof(double));
  //  for(i = 0; i < nCandMatch; i++) {
		//if((*newCluster).tempProb[i]==0)
		//	tranProb[i] = 0;
		//else
  //      tranProb[i] = exp( (*newCluster).trajSum[i] / getTemperT() );
  //      //tranProb[i] = exp( (*newCluster).trajSum[i] / getTemperT() ) * (double)((*newCluster).tempProb[i]);
  //      //tranProb[i] = ((*newCluster).sumAffinity + (*newCluster).trajSum[i])/getTemperT() * (double)((*newCluster).tempProb[i]);
  //  }
	memset(tranProb, 0, nCandMatch*sizeof(double));
	//for(i = 0; i < nWindow; i++) {
	int nWindowNext = nWindow * next;
	int nextnCandMatch = next * nCandMatch;
	for(i = 0; i < nWindow; i++) {
		int temp = rankM[i + nWindowNext];
		if((*newCluster).tempProb[temp]!=0) {
			tranProb[temp] = exp( M[nextnCandMatch + temp ]/getTemperT() );
		}
		//if((*newCluster).tempProb[rankM[i + nWindow*next]]==1) {
		//	tranProb[ rankM[i + nWindow*next] ] = exp(M[ rankM[i + nWindow*next] + nCandMatch*next ]/getTemperT());
		//}
	}

    for(i = 0; i < nCandMatch; i++) {
        sum += tranProb[i];
    }
    if(sum == 0) {
        for(i = 0; i < nCandMatch; i++) {
            tranProb[i] = (*newCluster).tempProb[i];
        }
    }
    if( !(*newCluster).termFlag ) {
        normalizeWeights( tranProb, nCandMatch );
    }
    
    // set the cumulative transition probability
    (*newCluster).cumTranProb = (double*)malloc(nCandMatch * sizeof(double));
    (*newCluster).cumTranProb[0] = tranProb[0];
    for(i = 1 ; i < nCandMatch; i++) {
        (*newCluster).cumTranProb[i] = (*newCluster).cumTranProb[i-1] + tranProb[i];
    }
    free(tranProb);
    
    double lastComponent = (*newCluster).cumTranProb[nCandMatch-1];
    for(j = nCandMatch-1; j >=0; j--) {
        if( lastComponent > (*newCluster).cumTranProb[j] ) {
            break;
        }
    }
    tempLastId = j;
    for(i = tempLastId+1; i< nCandMatch; i++) {
        (*newCluster).cumTranProb[i] = 1;
    }
    
    free(MNext);
}


/* set cluster */
void setCluster2(int next, int nMembers, cluster *newCluster, cluster *prevCluster, int *tempTraj) {

	//double *MNext = (double*)malloc(nCandMatch *sizeof(double));
    int i, j, tPoint;
    int tempProbSum = 0;
	double *tranProb;
    double sum = 0;
    int tempLastId;
	double lastComponent;

    // set the number of members
    (*newCluster).nMembers = nMembers;
    
    // set the trajectory
    (*newCluster).trajectory = (int*)malloc(l * sizeof(int));

	///*
	for(tPoint = 0; tPoint < l-1; tPoint++) {
		if(next < (*prevCluster).trajectory[tPoint]) break;
        else    (*newCluster).trajectory[tPoint] = (*prevCluster).trajectory[tPoint];
    }
    (*newCluster).trajectory[tPoint] = next;
    tPoint++;
    for(; tPoint < l; tPoint++) (*newCluster).trajectory[tPoint] = (*prevCluster).trajectory[tPoint-1];
    //*/
		
	int nextnCandMatch = next * nCandMatch;
    // set the sumAffinity
    //(*newCluster).sumAffinity = (*prevCluster).sumAffinity;
    //for(i = 0; i < l-1; i++) {
    //    (*newCluster).sumAffinity += M[nextnCandMatch + (*prevCluster).trajectory[i]];
    //}
	(*newCluster).sumAffinity = (*prevCluster).sumAffinity + (*prevCluster).trajSum[next];

    // set the transition probability
    (*newCluster).trajSum = (double*)malloc(nCandMatch * sizeof(double));
    for(i = 0; i < nCandMatch; i++) {
        //(*newCluster).trajSum[i] = (*prevCluster).trajSum[i] + MNext[i];
        (*newCluster).trajSum[i] = (*prevCluster).trajSum[i] + M[nextnCandMatch + i];
    }
	
    // set the tempProb
    (*newCluster).tempProb = (int*)malloc(nCandMatch * sizeof(int));
    for(i = 0; i < nCandMatch; i++) {
        (*newCluster).tempProb[i] = (*prevCluster).tempProb[i];
    }

    for(i = 0; i < confInfo.groupList1[confInfo.groupRef[next]].nMembers; i++) {
        (*newCluster).tempProb[ confInfo.groupList1[confInfo.groupRef[next]].members[i] ] = 0;
    }
    for(i = 0; i < confInfo.groupList2[confInfo.groupRef[next + nCandMatch]].nMembers; i++) {
        (*newCluster).tempProb[ confInfo.groupList2[confInfo.groupRef[next + nCandMatch]].members[i] ] = 0;
    }
    
    (*newCluster).termFlag = false;
    for(i = 0;i < nCandMatch; i++) {
        tempProbSum += (*newCluster).tempProb[i];
    }
    if(tempProbSum == 0)    (*newCluster).termFlag = true;    

    tranProb = (double*)malloc(nCandMatch * sizeof(double));
	memset(tranProb, 0, nCandMatch*sizeof(double));

	//for(i = 0; i < nWindow; i++) {
	//	int temp = rankM[i + nWindow*next];
	//	if((*newCluster).tempProb[temp]!=0) {
	//		for(j = 0; j < l-1; j++) {
	//			tranProb[temp] += M[ temp + nCandMatch*(*prevCluster).trajectory[j] ];
	//		}
	//		tranProb[ temp ] += M[nextnCandMatch + temp];
	//		tranProb[ temp ] = exp( tranProb[ temp ] / getTemperT() );
	//	}
	//}

	int nWindowNext = nWindow * next;
	for(i = 0; i < nWindow; i++) {
		int temp = rankM[i+nWindowNext];
		if((*newCluster).tempProb[temp]!=0)
		//tranProb[temp] = exp( (*newCluster).trajSum[temp] / getTemperT() );
		tranProb[temp] = kernel( (*newCluster).trajSum[temp] / getTemperT() );
    }

    for(i = 0; i < nCandMatch; i++) {
        sum += tranProb[i];
    }
    if(sum == 0) {
        for(i = 0; i < nCandMatch; i++) {
            tranProb[i] = (*newCluster).tempProb[i];
        }
    }


    normalizeWeights( tranProb, nCandMatch );

    // set the cumulative transition probability
    (*newCluster).cumTranProb = (double*)malloc(nCandMatch * sizeof(double));
    (*newCluster).cumTranProb[0] = tranProb[0];
    for(i = 1 ; i < nCandMatch; i++) {
        (*newCluster).cumTranProb[i] = (*newCluster).cumTranProb[i-1] + tranProb[i];
    }
    free(tranProb);
    
    lastComponent = (*newCluster).cumTranProb[nCandMatch-1];
    for(j = nCandMatch-1; j >=0; j--) {
        if( lastComponent > (*newCluster).cumTranProb[j] ) {
            break;
        }
    }
    tempLastId = j;
    for(i = tempLastId+1; i< nCandMatch; i++) {
        (*newCluster).cumTranProb[i] = 1;
    }
    
    //free(MNext);

}


void setGroupList(int *tempGroup1, int *tempGroup2) {
    int i, j;
    int mIndex;
	int nCandMatchI;

    for(i = 0; i < confInfo.nGroup1; i++) {
		nCandMatchI = nCandMatch * i;
        // get the number of members
        confInfo.groupList1[i].nMembers = 0;
        for(j = 0; j < nCandMatch; j++) {
            confInfo.groupList1[i].nMembers += (int)(tempGroup1[j + nCandMatchI]);
        }
        // save the group members
        confInfo.groupList1[i].members = (int*)malloc( confInfo.groupList1[i].nMembers * sizeof(int) );
        mIndex = 0;
        for(j = 0; j < nCandMatch; j++) {
            if( (int)(tempGroup1[j + nCandMatchI]) == 1) {
                confInfo.groupList1[i].members[mIndex] = j;
                mIndex++;
                confInfo.groupRef[j] = i;
            }
        }
    }
	
    for(i = 0; i < confInfo.nGroup2; i++) {
		nCandMatchI = nCandMatch * i;
        // get the number of members
        confInfo.groupList2[i].nMembers = 0;
        for(j = 0; j < nCandMatch; j++) {
            confInfo.groupList2[i].nMembers += (int)(tempGroup2[j + nCandMatchI]);
        }
        // save the group members
        confInfo.groupList2[i].members = (int*)malloc( confInfo.groupList2[i].nMembers * sizeof(int) );
        mIndex = 0;
        for(j = 0; j < nCandMatch; j++) {
            if( (int)(tempGroup2[j + nCandMatchI]) == 1) {
                confInfo.groupList2[i].members[mIndex] = j;
                mIndex++;
                confInfo.groupRef[j+nCandMatch] = i;
            }
        }
    }

	for(i = 0; i < 2*nCandMatch; i++) {
		if(confInfo.groupRef[i] == -1) {
			printf("error in groupRef\n");
			break;
		}
	}
}


/* free all components in clusters */
void freeClusters(cluster *clusters, int nCluster) {
    int i;
    for(i = 0; i < nCluster; i++) {
		if(clusters[i].trajectory != NULL)	free(clusters[i].trajectory);
        if(clusters[i].trajSum != NULL)		free(clusters[i].trajSum);
        if(clusters[i].tempProb != NULL)	free(clusters[i].tempProb);
        if(clusters[i].cumTranProb != NULL)	free(clusters[i].cumTranProb);
    }
}


void freeGroupList() {
    int i;
    for(i = 0; i < confInfo.nGroup1; i++) {
        free(confInfo.groupList1[i].members);
    }
    for(i = 0; i < confInfo.nGroup2; i++) {
        free(confInfo.groupList2[i].members);
    }
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


void desample(cluster* newClusters, cluster* newClustersDD) {
    int i, j, k, theMemberId, theMember;
    double minSum;
	double temp;

    for(i = 0; i < nNewCluster; i++) {
        
        // what member to remove?
        theMemberId = 0;
		theMember = newClusters[i].trajectory[0];
        minSum = newClusters[i].trajSum[ newClusters[i].trajectory[0] ];
        //minSum = 0;
		//for(k = 0; k < l; k++) minSum += M[ newClusters[i].trajectory[0] + nCandMatch * newClusters[i].trajectory[k] ];
        for(j = 1; j < l; j++) {
			//temp = 0;
			//for(k = 0; k < l; k++) temp += M[ newClusters[i].trajectory[j] + nCandMatch * newClusters[i].trajectory[k] ];
            temp = newClusters[i].trajSum[ newClusters[i].trajectory[j] ];
            if( temp < minSum ) {
                theMemberId = j;
                theMember = newClusters[i].trajectory[j];
                minSum = temp;
            }
        }
        
        desetCluster(theMember, theMemberId, minSum, i, newClusters, newClustersDD);
    }
}

void desetCluster(int theMember, int theMemberId, double minSum, int i, cluster *newClusters, cluster *newClustersDD) {
   
	// set MDel
    int ii, j, k, first;
	int theGroup1, theGroup2;
    double *MDel = (double*)malloc(nCandMatch *sizeof(double));
    bool tempE;
    int tempProbSum = 0;
	double *tranProb;
    double sum = 0;
    int tempLastId;
	double lastComponent;

	int theMembernCandMatch = theMember * nCandMatch;
	for(j = 0; j < nCandMatch; j++) {
        MDel[j] = M[theMembernCandMatch + j];
    }

    // update nMembers
    newClustersDD[i].nMembers = newClusters[i].nMembers;

    // update trajectory
    newClustersDD[i].trajectory = (int*)malloc( (l-1) * sizeof(int) );
    for(j = 0; j < theMemberId; j++) {
        newClustersDD[i].trajectory[j] = newClusters[i].trajectory[j];
    }
    for(; j < l-1; j++) {
        newClustersDD[i].trajectory[j] = newClusters[i].trajectory[j+1];
    }

    // update sumAffinity
    newClustersDD[i].sumAffinity = newClusters[i].sumAffinity - minSum;

    // set the transition probability
    newClustersDD[i].trajSum = (double*)malloc(nCandMatch * sizeof(double));
    for(j = 0; j < nCandMatch; j++) {
        newClustersDD[i].trajSum[j] = newClusters[i].trajSum[j] - MDel[j];
    }

    // update the tempProb
    newClustersDD[i].tempProb = (int*)malloc(nCandMatch * sizeof(int));
    for(j = 0; j < nCandMatch; j++) {
        newClustersDD[i].tempProb[j] = newClusters[i].tempProb[j];
    }
    
    theGroup1 = confInfo.groupRef[theMember];
    theGroup2 = confInfo.groupRef[theMember + nCandMatch];
    for(j = 0; j < confInfo.groupList1[theGroup1].nMembers; j++) {
        tempE = false;
        for(k = 0; k < l-1; k++) {
            if( confInfo.groupRef[ newClustersDD[i].trajectory[k] + nCandMatch] == confInfo.groupRef[ confInfo.groupList1[theGroup1].members[j] + nCandMatch ] ) {
                tempE = true;
                break;
            }
        }
        if(tempE == false) {
            newClustersDD[i].tempProb[ confInfo.groupList1[theGroup1].members[j] ] = 1;
        }
    }
    for(j = 0; j < confInfo.groupList2[theGroup2].nMembers; j++) {
        tempE = false;
        for(k = 0; k < l-1; k++) {
            if( confInfo.groupRef[ newClustersDD[i].trajectory[k] ] == confInfo.groupRef[ confInfo.groupList2[theGroup2].members[j] ] ) {
                tempE = true;
                break;
            }
        }
        if(tempE == false) {
            newClustersDD[i].tempProb[ confInfo.groupList2[theGroup2].members[j] ] = 1;
        }
    }

    // update termFlag
    newClustersDD[i].termFlag = false;
    for(j = 0; j < nCandMatch; j++) {
        tempProbSum += newClustersDD[i].tempProb[j];
    }
    if(tempProbSum == 0) {
        newClustersDD[i].termFlag = true;
        printf("something wrong. in the desample function\n");
    }

    // set tranProb to calculate cumTranProb
    tranProb = (double*)malloc(nCandMatch * sizeof(double));
	memset(tranProb, 0, nCandMatch*sizeof(double));
	first = newClustersDD[i].trajectory[0];
	int nWindowFirst = nWindow*first;
    for(j = 0; j < nWindow; j++) {
		int temp = rankM[j + nWindowFirst];
		if((newClustersDD[i]).tempProb[temp]==0)
			tranProb[temp] = 0;
		else
			//tranProb[temp] = exp( newClustersDD[i].trajSum[temp] / getTemperT() );
			tranProb[temp] = kernel( newClustersDD[i].trajSum[temp] / getTemperT() );
    }
	//memset(tranProb, 0, nCandMatch*sizeof(double));
	//first = newClustersDD[i].trajectory[0];
	//for(ii = 0; ii < nWindow; ii++) {
	//	int temp = rankM[ii + nWindow*first];
	//	if(newClustersDD[i].tempProb[temp]!=0) {
	//		for(j = 0; j < l-1; j++) {
	//			tranProb[ temp ] += M[ temp + nCandMatch*newClustersDD[i].trajectory[j] ];
	//		}
	//		tranProb[ temp ] = exp(tranProb[ temp ]/getTemperT());
	//	}
	//}

	for(j = 0; j < nCandMatch; j++) {
        sum += tranProb[j];
    }
    if(sum == 0) {
        for(j = 0; j < nCandMatch; j++) {
            tranProb[j] = newClustersDD[i].tempProb[j];
        }
    }
    if( !newClustersDD[i].termFlag ) {
        normalizeWeights( tranProb, nCandMatch );
    }

    // update the cumulative transition probability
    newClustersDD[i].cumTranProb = (double*)malloc(nCandMatch * sizeof(double));
    newClustersDD[i].cumTranProb[0] = tranProb[0];
    for(j = 1 ; j < nCandMatch; j++) {
        newClustersDD[i].cumTranProb[j] = newClustersDD[i].cumTranProb[j-1] + tranProb[j];
    }
    free(tranProb);

    lastComponent = newClustersDD[i].cumTranProb[nCandMatch-1];
    for(j = nCandMatch-1; j >=0; j--) {
        if( lastComponent > newClustersDD[i].cumTranProb[j] ) {
            break;
        }
    }
    tempLastId = j;
    for(k = tempLastId+1; k< nCandMatch; k++) {
        newClustersDD[i].cumTranProb[k] = 1;
    }


    free(MDel);
}


// return how many particles exit
//int makeClusterId(int* clusterId, cluster* newClusters, cluster** saveCluster, double* nTrajectory, double* trajectoryIdx, int* trajCount, int* nTrajCount) {
int makeClusterId(int* clusterId, cluster* newClusters, cluster** saveCluster, double* nTrajectory, double* trajectoryIdx, int* trajCount, int* nTrajCount) {

	int i, j, k(0), nExit(0);
    
    for(i = 0; i < nNewCluster; i++) {
        if(newClusters[i].termFlag) {
            nExit += newClusters[i].nMembers;
            *(nTrajectory+(*nTrajCount)) = l;
            for(j = 0; j < l; j++) {
                *(trajectoryIdx+(*trajCount)+j) = newClusters[i].trajectory[j];
            }
            (*nTrajCount)++;
            (*trajCount) += l;
            
            //save
            if(newClusters[i].sumAffinity > (**saveCluster).sumAffinity) {
                freeClusters((*saveCluster),1);
                (*saveCluster) = (cluster*)malloc(sizeof(cluster));
                (**saveCluster).sumAffinity = newClusters[i].sumAffinity;
                (**saveCluster).trajectory = (int*)malloc(sizeof(int)*l);
                for(j = 0; j < l; j++) {
                    (**saveCluster).trajectory[j] = newClusters[i].trajectory[j];
                }
                (**saveCluster).trajSum = NULL;
                (**saveCluster).tempProb = NULL;
                (**saveCluster).cumTranProb = NULL;
				(**saveCluster).length = l;
            }
        }
        else {
            for(j = 0; j < newClusters[i].nMembers; j++) {
                clusterId[k] = i;
                k++;
            }
        }
    }
    return nExit;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////



/* sample one particle from cumProb */
int sample(double* cumProb, int n) {
//     if(1-cumProb[n-1]>0.00001)   {
//         printf("invalid cum_prob\n");
//     }
	int i;
     //rand()
     double u = (double)rand() / (double)RAND_MAX;
    //double u = Random();
    if(u == 0) {
        for(i = 0; i < n; i++) {
            if(u < cumProb[i]) return i;
        }
    }
    else {
        for(i = 0; i < n; i++) {
            if(u <= cumProb[i]) return i;
        }
    }
    return (n-1);
    
//     // Random()
//     double u = Random();
//     for(i = 0; i < n; i++) {
//             if(u <= cumProb[i]) return i;
//     }
//     return (n-1);
    
//     printf("you should not be in here.\n");
}



/* normalize x (make sum(x) = 1) */
void normalizeWeights(double *x, int nX) {

	int i;
	double sum;

	if(nX==1) {
		x[0] = 1;
	}
	else{
		sum = 0;
		for(i = 0; i < nX; i++) {
			sum += x[i];
		}
		if(sum == 0) {
// 	        printf("you should not be in here\n");
			return;
		}
		for(i = 0; i < nX; i++) {
			x[i] /= sum;
		}
	}

}


/* count the number of unique elements */
int countUnique(int *x, int nX, int **nMembers) {   
    int i; 
    int nUnique = 1;
    int *tempNMembers = (int*)malloc( nX * sizeof(int) );
	memset(tempNMembers, 0, nX * sizeof(int));

// x is already sorted
//     qsort(x);
    tempNMembers[0] = 1;
    for(i = 1; i < nX; i++) {
        if(x[i] != x[i-1])  nUnique++;
        tempNMembers[nUnique-1]++;
    }
    
    // output assign
    *nMembers = (int*)malloc( nUnique * sizeof(int));
    for(i = 0; i < nUnique; i++){
        *(*nMembers + i) = tempNMembers[i];
    }
    free(tempNMembers);
    return nUnique;
}


/* compare two numbers (for qsort) */
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}



bool isEqualComb(int *comb1, int *comb2) {
    int i = 0;
    for( i = 0; i < l; i++) {
        if(comb1[i] != comb2[i])    return false;
    }
    return true;
}


int* radixSort(int* tempTraj) {
	int i, tt;
	int count(0);
	LLNode *buckets, *tempBuckets, **bucketEnds, **tempBucketEnds;
	buckets = (LLNode*)malloc(nCandMatch * sizeof(LLNode));
	for(i = 0; i < nCandMatch; i++)	buckets[i].next = NULL;
	tempBuckets = (LLNode*)malloc(nCandMatch * sizeof(LLNode));
	for(i = 0; i < nCandMatch; i++)	tempBuckets[i].next = NULL;
	bucketEnds = (LLNode**)malloc(nCandMatch * sizeof(LLNode*));
	for(i = 0; i < nCandMatch; i++)	bucketEnds[i] = buckets+i;
	tempBucketEnds = (LLNode**)malloc(nCandMatch * sizeof(LLNode*));
	for(i = 0; i < nCandMatch; i++)	tempBucketEnds[i] = tempBuckets+i;

	int bucketIdx;
	int tempParticleId;
	LLNode botton;
	int *sortParticleId;

	// sort particles
	for(i = 0; i < nParticle; i++) {
		bucketIdx = tempTraj[i];
		LLInsert(bucketEnds, bucketIdx, i);
	}

	for(tt = 1; tt < l; tt++) {
		// initialize
		tempBuckets = (LLNode*)malloc(nCandMatch *sizeof(LLNode));
		for(i = 0; i < nCandMatch; i++) tempBuckets[i].next = NULL;
		for(i = 0; i < nCandMatch; i++) tempBucketEnds[i] = tempBuckets+i;
				
		// distribute
		for(i = 0; i < nCandMatch; i++) {
			if(buckets[i].next != NULL) {
				botton = *(buckets[i].next);
				while(1){
					tempParticleId = botton.data;
					bucketIdx = tempTraj[tempParticleId + tt*nParticle];
					
					LLInsert(tempBucketEnds, bucketIdx, tempParticleId);
					
					if(botton.next == NULL)	break;
					else botton = *(botton.next);
				}
			}
		}

		// update buckets
		freeLL(buckets);
		buckets = tempBuckets;
		tempBuckets = NULL;
	}

	count = 0;
	sortParticleId = (int*)malloc(nParticle * sizeof(int));
    for(i = 0; i < nCandMatch; i++) {
		botton = buckets[i];
		while(botton.next!=NULL) {
			botton = *(botton.next);
            sortParticleId[count] = botton.data;
			count++;
		}
    }
	freeLL(buckets);
	free(bucketEnds);
	free(tempBucketEnds);
	return sortParticleId;	// TODO: free sortParticleId
}

void LLInsert(LLNode **bucketEnds, int bucketIdx, int data) {
	LLNode *newNode;
	newNode = (LLNode*)malloc(sizeof(LLNode));
	newNode->data = data;
	newNode->next = NULL;

	bucketEnds[bucketIdx]->next = newNode;
	bucketEnds[bucketIdx] = newNode;
	newNode = NULL;
}

void freeLL(LLNode *buckets) {
	int i;
	LLNode *delNode, *temp;
	
	for(i = 0; i < nCandMatch; i++) {
		temp = (buckets+i)->next;
		while(temp != NULL) {
			delNode = temp;
			temp = temp->next;
			free(delNode);
		}
	}
	free(buckets);
}


void insertNew(int xPfPred, int tempClusterId, int* tempTraj, int j) {
    int tPoint;
    for(tPoint = 0; tPoint < l-1; tPoint++) {
        if(xPfPred < clusters[tempClusterId].trajectory[tPoint]) break;
        else    tempTraj[j + tPoint*nParticle] = clusters[tempClusterId].trajectory[tPoint];
    }
    tempTraj[j + tPoint*nParticle] = xPfPred;
    tPoint++;
    for(; tPoint < l; tPoint++) {
		tempTraj[j + tPoint*nParticle] = clusters[tempClusterId].trajectory[tPoint-1];
	}
}

double kernel(double x) {
	double y;
	y = exp(x);
	return y;
}

double getTemperT() {
    //double temperT;
    // temperT = 2;
    //return temperT;
	return temperature;
}