 #include "particles_eccv.h"

// Solve graph matching problem using Sequential Monte Carlo
//
// INPUT :  Affinity matrix		(nCandMatch x nCandMatch)
//			confliction group 1 (nCandMatch x nGroup1)
//			confliction group 2 (nCandMatch x nGroup2)
//			rankM				(unused)
//       # of particles		(scalar)
//			tau					(scalar)
// OUTPUT : solution (nCandMatch x 1) binary indicator vector

confGroupInfo confInfo;
int nCandMatch;			// number of candidate matches
double *M;				// affinity matrix
int nParticle;			// number of particles
int l;					// number of matches in each particle
cluster *clusters;
int nCluster, nNewCluster;
int *rankM;				// unused
int nWindow;			// unused
double temperature;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double MAX_DOUBLE = 1.7e308;
    const double MIN_DOUBLE = 1.7e-308;
    const int MAX_L = 1000;
    const int MAX_MATCHES = 500000;
	
	int i, j, k, t(1);		// counting variables
    mxClassID tempCategory;	// type variable

	// number variables
	int nExit(0);
	int nNewGroup;
	int *nGroupMembers;

	// particle variables
    int *M_init, *M_samples, *M_particles;		// members are indices wrt match
    int *P_rspl, *P_sortPtcls;					// members are indices wrt particle
	int *C_particles, *C_prevPtcls, C_temp;		// members are indices wrt cluster

	// cluster variables
	cluster *newClusters, *newClustersDD, *saveCluster;
    int *MC_news;
    int *groupStartPoints;
	double *groupWeights, *cumGroupWeights;
	//
	double tempSumAffinity;
	double tempSum_;

    bool terminate = false;

    /* Input */
    //double *M;		// affinity matrix
    int *nParticlePr;
    int *tempGroup1;
	int *tempGroup2;
	double *temperaturePr;

	// TODO: check input type
    M = mxGetPr(prhs[0]);
    tempGroup1 = (int32_T *)mxGetData(prhs[1]);
    tempGroup2 = (int32_T *)mxGetData(prhs[2]);
	rankM = (int32_T *)mxGetData(prhs[3]); // unused
	nParticlePr = (int32_T *)mxGetData(prhs[4]);
	temperaturePr = mxGetPr(prhs[5]);

	nParticle = *nParticlePr;
    nCandMatch = mxGetM(prhs[0]);
    confInfo.nGroup1 = mxGetN(prhs[1]);
    confInfo.nGroup2 = mxGetN(prhs[2]);
	nWindow = mxGetM(prhs[3]);		// unused
	temperature = *temperaturePr;

	confInfo.groupList1 = (groupList*)malloc(confInfo.nGroup1 * sizeof(groupList));
    confInfo.groupList2 = (groupList*)malloc(confInfo.nGroup2 * sizeof(groupList));
    confInfo.groupRef = (int*)malloc(2 * nCandMatch * sizeof(int));
    memset(confInfo.groupRef, -1, 2*nCandMatch*sizeof(int));
	setGroupList(tempGroup1, tempGroup2);
	
	/* Output */
    double *X;	// TODO: change type from double to int32
    plhs[0] = mxCreateDoubleMatrix(nCandMatch, 1,mxREAL);
    X = mxGetPr(plhs[0]);
    
    // @
    double* nTrajectory;
	double* trajectoryIdx;
    int trajCount = 0;
    int nTrajCount = 0;
    plhs[1] = mxCreateDoubleMatrix(10000, 1, mxREAL);
    nTrajectory = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(MAX_MATCHES, 1, mxREAL);
    trajectoryIdx = mxGetPr(plhs[2]);

	srand(time(NULL));
	//srand(0);
    
    ///////////////////////////////////////////////////////////////////////
    /* Initialization */
	t = 1;
	l = 1;
	nCluster = 1;
    initDistribution(&M_init);
    
    // initialize cluster id for each particle
    C_particles = (int*)malloc(nParticle * sizeof(int));
	memset(C_particles, 0, nParticle*sizeof(int));
        
	saveCluster = (cluster*)malloc(sizeof(cluster));
    (*saveCluster).sumAffinity = -1;
    (*saveCluster).trajectory = NULL;
    (*saveCluster).trajSum = NULL;
    (*saveCluster).tempProb = NULL;
    (*saveCluster).cumTranProb = NULL;
    
    group(M_init, &newClusters);
    nExit = makeClusterId(C_particles, newClusters, &saveCluster, nTrajectory, trajectoryIdx, &trajCount, &nTrajCount);
	if(nExit != 0) {
        printf("somethig wrong!\n");
    }
    
    clusters = newClusters;
    newClusters = NULL;
    nCluster = nNewCluster;
    free(M_init);

    
    ///////////////////////////////////////////////////////////////////////
    /* SMC loop */
    M_samples = (int*)malloc(nParticle * sizeof(int));
	while(1) {
        
        t++;
        
        /* termination condition */
        if(nParticle <= 0) {
            break;
        }
        if( l > MAX_L ) {
            printf("something wrong. %d-th iteration \n",l);
            break;
        }
        l++;
        
        /* perform prediction and measurement for each particle */        
        M_particles = (int*)malloc( l * nParticle * sizeof(int) );   // *
        for(j = 0; j < nParticle; j++) {
            // 1. Sample x
            C_temp = C_particles[j];
            M_samples[j] = sample(clusters[C_temp].cumTranProb, nCandMatch);
            insertNew(M_samples[j], C_temp, M_particles, j);
        }
		
		// sort particles
        P_sortPtcls = radixSort(M_particles);
        
        nNewGroup = 1;
        groupStartPoints = (int*)malloc(nParticle * sizeof(int));
        C_prevPtcls = (int*)malloc(nParticle * sizeof(int));
        MC_news = (int*)malloc(nParticle * sizeof(int));
        groupStartPoints[0] = 0;
		C_prevPtcls[0] = C_particles[P_sortPtcls[0]];
		MC_news[0] = M_samples[ P_sortPtcls[0] ];
        for(i = 1; i < nParticle; i++) {
            int *M_tmpPtcl1 = (int*)malloc(l * sizeof(int));
            int *M_tmpPtcl2 = (int*)malloc(l * sizeof(int));
            for(j = 0; j < l; j++) {
                M_tmpPtcl1[j] = M_particles[P_sortPtcls[i-1] + j*nParticle];
                M_tmpPtcl2[j] = M_particles[P_sortPtcls[i] + j*nParticle];
            }
            
            if( !isEqualComb(M_tmpPtcl1, M_tmpPtcl2) ) {
                C_prevPtcls[nNewGroup] = C_particles[P_sortPtcls[i]];
                MC_news[nNewGroup] = M_samples[ P_sortPtcls[i] ];
                groupStartPoints[nNewGroup] = i;
                nNewGroup++;
            }
            free(M_tmpPtcl1);
            free(M_tmpPtcl2);
        }
        
        // update nGroupMembers
        nGroupMembers = (int*)malloc(nNewGroup * sizeof(double));
        for(i = 0; i < nNewGroup-1; i++) {
            nGroupMembers[i] = groupStartPoints[i+1] - groupStartPoints[i];
        }
        nGroupMembers[nNewGroup-1] = nParticle - groupStartPoints[i];
        
        // update groupWeights
        groupWeights = (double*)malloc(nNewGroup * sizeof(double));

		double avgSum = 0;
		for(i = 0; i < nNewGroup; i++) {
			avgSum += clusters[C_prevPtcls[i]].sumAffinity;
		}
		avgSum = avgSum / nNewGroup;

        for(i = 0; i < nNewGroup; i++) {
            tempSumAffinity = clusters[C_prevPtcls[i]].sumAffinity - avgSum;

			for(j = 0; j < l-1; j++) {
                tempSumAffinity += M[MC_news[i] + nCandMatch*(clusters[C_prevPtcls[i]].trajectory[j])];
            }
            
            if( (tempSumAffinity / getTemperT()) > log(MAX_DOUBLE / nNewGroup) ) {
				printf("numerical problem1 in SMCM-new\n");
                groupWeights[i] = MAX_DOUBLE / nNewGroup;
            }
            else if( (tempSumAffinity / getTemperT()) < log(MIN_DOUBLE) ) {
                groupWeights[i] = 0;
            }
            else {
                groupWeights[i] = exp(tempSumAffinity/getTemperT());
            }
        }

        // @ adhoc
            tempSum_ = 0;
            for( i = 0; i < nNewGroup; i++) {
                tempSum_ += groupWeights[i];
            }
            if( tempSum_ < 0.00001 ) {
				printf("numerical problem2 in SMCM-new\n");
                for(i = 0; i < nNewGroup; i++) {
                    groupWeights[i] = 1/nNewGroup;
                }
            }
            else {
                // normalize
                normalizeWeights(groupWeights, nNewGroup);
            }

        
        // get cumulative weight distribution
        cumGroupWeights = (double*)malloc(nNewGroup * sizeof(double));
        cumGroupWeights[0] = groupWeights[0];
        for(i = 1; i < nNewGroup; i++) {
            cumGroupWeights[i] = cumGroupWeights[i-1] + groupWeights[i];
        }
                
        // 2. Resample
        P_rspl = (int*)malloc(nParticle * sizeof(int));
        for(i = 0; i < nParticle; i++) {
            P_rspl[i] = sample(cumGroupWeights,nNewGroup);
        }

        group2(P_rspl, &newClusters, C_prevPtcls, MC_news, M_particles, P_sortPtcls);

        // make clusterId
        nExit = makeClusterId(C_particles, newClusters, &saveCluster, nTrajectory, trajectoryIdx, &trajCount, &nTrajCount);
		nParticle = nParticle - nExit;
		
		// free memory
        free(P_sortPtcls);
        free(P_rspl);
        free(M_particles);
        free(MC_news);
        free(C_prevPtcls);
        free(groupStartPoints);
        free(nGroupMembers);
        free(groupWeights);
        free(cumGroupWeights);

		// 3. Update
		if( t%2 ==1 ) {
           newClustersDD = (cluster*)malloc(nNewCluster * sizeof(cluster));
           desample(newClusters, newClustersDD);
           l--;

           freeClusters(clusters, nCluster);
           free(clusters);
           freeClusters(newClusters, nNewCluster);
           free(newClusters);
           clusters = newClustersDD;
           newClustersDD = NULL;
           nCluster = nNewCluster;
       }
       else {

            freeClusters(clusters, nCluster);
            free(clusters);
            clusters = newClusters;
            newClusters = NULL;
            nCluster = nNewCluster;
			//printf("nCluster = %d\n", nCluster);
			//printf("nNewGroup = %d\n", nNewGroup);
	        }
    }
        
    /* output assigning */
	memset(X, 0, nCandMatch*sizeof(double));
    for(i = 0; i < (*saveCluster).length; i++) {
        X[ (*saveCluster).trajectory[i] ] += 1;
    }

   
    free(M_samples);
    free(C_particles);
	freeGroupList();
}