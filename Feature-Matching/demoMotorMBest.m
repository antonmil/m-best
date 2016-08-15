%% Code for comparison between different graph matching solvers (QBP program)
%  on Pascal VOC matching datasetes (MotorBike dataset) and their improved
%  results using m-best solutions.  

% The code is tested under Windows 10 (64bit)

% If you use this code for your research, please, cite:
% S. H. Rezatofighi, A. Milan, Z. Zhang, Q. Shi, A. Dick, I. Reid, "Joint 
% Probabilistic Matching Using m-Best Solutions ", IEEE Conference 
% on Computer Vision and Pattern Recognition (CVPR), 2016.

% Note: The results may be slightly different from the reported results due
% to some randomnesses.

% Author: S. Hamid Rezatofighi
% Last updated: July 30, 2016

% For questions contact the author at: hamid.rezatofighi@adelaide.edu.au

clear variables;
global footpath;
footpath = cd;

%% Experiment setup
rng(45678);
AddPath % Add the necessary paths
prSet(1);

MaxInstance = 20; % Number of instance (for Car= 30, for Motorbike = 20)
NofAlgos = 11; % Number of algorithms
MaxOutliers = 20; % Maximum number of outliers added to the problem
parKnl = st('alg', 'pas'); % type of affinity: only edge distance

%% src parameter
tag = 'pas';

%% BP parameters
bpoptions.outIter = 1;
bpoptions.innerIter = 10;
BaBoptions.MaxIter = 100;
BaBoptions.bpoptions = bpoptions;

%% Mbest parameters
param.mbest=5; % Number of solutions
param.chck_sols=1; 

%% Pre-allocations
AvgAcc = zeros(MaxOutliers+1, NofAlgos);
AvgObj = zeros(MaxOutliers+1, NofAlgos);
StdAcc = zeros(MaxOutliers+1, NofAlgos);
StdObj = zeros(MaxOutliers+1, NofAlgos);
AvgTime = zeros(MaxOutliers+1, NofAlgos);
StdTime = zeros(MaxOutliers+1, NofAlgos);
accs = cell(MaxOutliers+1, MaxInstance);
objs = cell(MaxOutliers+1, MaxInstance);
times = cell(MaxOutliers+1, MaxInstance);
sols = cell(MaxOutliers+1, MaxInstance);
m_sols = cell(MaxOutliers+1, MaxInstance);
m_objs = cell(MaxOutliers+1, MaxInstance);

for nOut = 0:MaxOutliers
    for kFs=1:MaxInstance
        disp(['Example #',num2str(kFs),' /Number of added outlier = ',num2str(nOut)])
        acc = zeros(1, NofAlgos);
        obj = zeros(1, NofAlgos);
        tm = zeros(1, NofAlgos);
        
        parKnl = st('alg', 'pas2'); % type of affinity: only edge distance
        %% algorithm parameter
        [pars, algs] = gmPar(2);
        par_mb=pars{6};par_mb{1,3}.alg='ipfp_mbst';
        par_mb{1,3}.mbst=param.mbest;
        par_mb{1,3}.chck_sols=param.chck_sols;
        %% src
        wsSrc = motorAsgSrc(kFs, nOut, 0);
        asgT = wsSrc.asgT;
        
        parG = st('link', 'del'); 
        parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
        wsFeat = motorAsgFeat(wsSrc, parG, parF, 'svL', 1);
        [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');
        
        [KP, KQ] = conKnlGphPQD(gphs, parKnl);
        K = conKnlGphKD(KP, KQ, gphs);
        Ct = ones(size(KP));
        
        %% SMCM (Sequential Monte Carlo Method)
        % problem = makeProblem_SMC( K,size(KP,1),size(KP,2),asgT.X);
        % setMethods;
        % [acc(1),obj(1),tm(1)] = wrapper_GM(methods(1), problem);
        
        %% Hypergraph
        % The results for Hypergraph method are extracted from the paper
        
        % Nguyen et al., "A flexible tensor block coordinate ascent scheme 
        % for hypergraph matching", CVPR 2015.
        %% PM
        tPm = tic;
        asgPm = pm(K, KQ, gphs, asgT);
        tm(2) = toc(tPm);
        acc(2) = asgPm.acc;
        obj(2) = asgPm.obj;
        sols{nOut+1,kFs}.PM=asgPm.X(:);

        %% GA
        tGa = tic;
        asgGa = gm(K, Ct, asgT, pars{1}{:});
        tm(3) = toc(tGa);
        acc(3) = asgGa.acc;
        obj(3) = asgGa.obj;
        sols{nOut+1,kFs}.GA=asgGa.X(:);
             
        %% SM
        tSm = tic;
        asgSm = gm(K, Ct, asgT, pars{3}{:});
        tm(4) = toc(tSm);
        acc(4) = asgSm.acc;
        obj(4) = asgSm.obj;
        sols{nOut+1,kFs}.SM=asgSm.X(:);
        
        
        %% SMAC
        tSmac = tic;
        asgSmac = gm(K, Ct, asgT, pars{4}{:});
        tm(5) = toc(tSmac);
        acc(5) = asgSmac.acc;
        obj(5) = asgSmac.obj;
        sols{nOut+1,kFs}.SMAC=asgSmac.X(:);
             
        %% IPFP-S
        tIpfp = tic;
        asgIpfpS = gm(K, Ct, asgT, pars{6}{:});
        tm(6) = toc(tIpfp);
        acc(6) = asgIpfpS.acc;
        obj(6) = asgIpfpS.obj;
        sols{nOut+1,kFs}.IPFPS=asgIpfpS.X(:);
        
        %% RRWM
        tRrwm = tic;
        asgRrwm = gm(K, Ct, asgT, pars{7}{:});
        tm(7) = toc(tRrwm);
        acc(7) = asgRrwm.acc;
        obj(7) = asgRrwm.obj;
        sols{nOut+1,kFs}.RRWM=asgRrwm.X(:);
        
        %% FGM-D
        tFgmD = tic;
        asgFgmD = fgmD(KP, KQ, Ct, gphs, asgT, pars{9}{:});
        tm(8) = toc(tFgmD);
        acc(8) = asgFgmD.acc;
        obj(8) = asgFgmD.obj;
        sols{nOut+1,kFs}.FGMD=asgFgmD.X(:);  
        %% BP MAP solver
        tbp = tic;
        asgBP = Hungarian_BP(K, Ct, asgT,BaBoptions);
        tm(9) = toc(tbp);        
        acc(9) = asgBP.acc;
        obj(9) = asgBP.obj;
        sols{nOut+1,kFs}.BP=asgBP.X(:); 
        
        %% M-Best IPFP-S
        tIpfp = tic;
        asgIpfpSMbst = gm(K, Ct, asgT, par_mb{:});
        tm(10) = toc(tIpfp);
        acc(10) = asgIpfpSMbst.acc;
        obj(10) = 0; % objective using the m-best approach (marginals) is not
        % comparable with the solvers objectives (estimated using MAP). 
        
        % storing all m solutions and their objectives
         sols{nOut+1,kFs}.IPFPMBST=asgIpfpSMbst.Xmbst(:);
         m_sols{nOut+1,kFs}.IPFPMBST=asgIpfpSMbst.X;
         m_objs{nOut+1,kFs}.IPFPMBST=asgIpfpSMbst.obj;

       %% M-Best BP 
        tbp = tic;
        asgBPMbst = mbest_BP(K, Ct, asgT,BaBoptions, param);
        tm(11) = toc(tbp);        
        acc(11) = asgBPMbst.acc;
        obj(11) = 0; % objective using the m-best approach (marginals) is not
        % comparable with the solvers objectives (estimated using MAP). 
        
        % storing all m solutions and their objectives
        sols{nOut+1,kFs}.BPMBST=asgBPMbst.Xmbst(:);
        m_sols{nOut+1,kFs}.BPMBST=asgBPMbst.X;
        m_objs{nOut+1,kFs}.BPMBST=asgBPMbst.obj;
        
        % print information
        times{nOut+1,kFs} = tm;
        accs{nOut+1,kFs} = acc;
        objs{nOut+1,kFs} = obj;
        save('Results\MotorResult_Mbst.mat');
      
    end
    tiimes = cell2mat(times(nOut+1,:)');
    objjs = cell2mat(objs(nOut+1,:)');
    acccs = cell2mat(accs(nOut+1,:)');
    maxobjs = max(objjs,[], 2);
    normalised_objs = zeros(MaxInstance,NofAlgos);
    for i=1:NofAlgos
        normalised_objs(:,i) = objjs(1, i)./maxobjs(1);
    end
    
    for i=1:NofAlgos
        AvgAcc(nOut+1, i) = mean(acccs(:,i));
        StdAcc(nOut+1, i) = std(acccs(:,i));
        AvgObj(nOut+1, i) = mean(normalised_objs(:, i));
        StdObj(nOut+1, i) = std(normalised_objs(:,i));
        AvgTime(nOut+1, i) = mean(tiimes(:,i));
        StdTime(nOut+1, i) = std(tiimes(:,i));
    end    
    save(['Results',filesep,'MotorResult_Mbst.mat']);
end

%% Plot Matching Accuracy
algs2 = algs([2 1 3:4,6:7 9]);
Algs = ['SMCM',algs2,'BP','mbst-IPFP','mbst-BP'];
stp=1;
colors = distinguishable_colors(length(Algs)-2);
figure;
for ji=1:length(Algs)
     if strcmp(Algs(ji),'mbst-BP')
        linstyle='--';
        clrr= colors(ji-2,:);
     elseif strcmp(Algs(ji),'mbst-IPFP')
          linstyle='--';
          clrr = colors(ji-4,:);
     else
         linstyle='-';
         clrr= colors(ji,:);
     end
    hold on, plot(0:stp:MaxOutliers, AvgAcc(1:stp:MaxOutliers+1,ji),linstyle, 'Color', clrr, 'LineWidth', 3)
end


set(gca,'FontName','Times','FontSize',24),
xlabel('#Outliers');
ylabel('Accuracy');
title('MotorBike') 

xlim([0,20]);
ylim([0.20,1]);
set(gca,'FontName','Times','FontSize',20),
h=legend(Algs,'Location','southwest');
set(gca,'FontName','Times','FontSize',20),
