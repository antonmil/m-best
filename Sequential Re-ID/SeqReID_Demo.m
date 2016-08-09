%% Code for SeqReID
 
% In this problem, the aim is to successfully track visually similar objects 
% in a video sequence by considering their appearance only and without
% using a motion prior. 
 
% This code uses the code package, experimental setup and the results published 
% by [Dicle et al. ICCV 2013]
 
% C. Dicle, O. I. Camps, and M. Sznaier. The way they move: Tracking 
% multiple targets with similar appearance. In ICCV, 2013.
 
% We have modified it by adding our SeqReID into this code package. 
% Similar to [Dicle et al. ICCV 2013], we use the annotated bounding box 
% provided for every target.

% Note 1: The results may be slightly different from the reported results due
% to some randomnesses.

 
% Note 2: The SeqReID is very simple tracker (Sequential matching only). It 
% is very robust to distinguish visually similar objects. However, if you  
% want to use it as a tracker working on real detections you may need to make 
% the tracker more robust or to incorporate the m-best idea in combination 
% with any other tracking algorithm.  This is mainly due to this fact that
% this simple tracker can be good to solve data association task and cannot 
% potentially deal with other multi-target tracking difficulties such as 
% false detections and target initiation & termination. 
 
 
% The code is tested under Windows 10 (64bit)
 
% If you use this code for your research, please, cite:
% S. H. Rezatofighi, A. Milan, Z. Zhang, Q. Shi, A. Dick, I. Reid, "Joint 
% Probabilistic Matching Using m-Best Solutions ", IEEE Conference 
% on Computer Vision and Pattern Recognition (CVPR), 2016.
 
 
% Author: S. Hamid Rezatofighi
% Last updated: August 5, 2016
 
% For questions contact the author at: hamid.rezatofighi@adelaide.edu.au


clear,clc,close all
AddPath % Add the necessary paths
rng(1233245);
%% Experiment setup

% Datasets to be tested
datasetPath = 'Data'; % dataset directory path
dataSets = {'crowd','seagulls','tud-crossing','tud-campus','PETS',...
    'AFL_S1','AFL_S2','AFL_S3'};

%% Parameters for smot tracker [Dicle et al. ICCV 2013] 
methods = {'ihtls','admm','ip'};
method  = methods{1};
saveOutput = false;

% Noise type to be tested
noiseType = {'fn','fp'};
noise_tag = noiseType{1};

% Noise levels
fn = [0:0];     
fp = [0.0:0.1:0.5];
N = length(fn);

% Number of trials per noise level.
TRIALS = 1;
% noise base
noisebase.fn = 0.0;
noisebase.fp = 0.0;
noisebase.gn = 0;
% mota base
motbase.fn = 0;
motbase.fp = 0;
motbase.mme = 0;
motbase.g = 0;
motbase.mota = 0;
motbase.mmerat = 0;

%% SeqReID parameters
param.min_dist=15; % Distance to missed detection
param.term_f=1; % Termination parameter (Terminate a track after 'term_f' 
% frames of missed detection) 
param.G=350;% Gate: Mamximum distance to consider matching candidates (pixels)
% This parameter is set to ease the computational complexity. Its lower 
% value makes the computaional burden significantly less and may improve 
% the results. But it would be against our "without a motion prior" assumption. 

%% Parameters for visualization
param.vis = false; % Visualsing trajectories

mot = cell(length(dataSets),2);
for ii = 1:length(dataSets)
    
    seqName = dataSets{ii}; % Dataset name
    LoadDetections; % Load the annotated detection boxes (idl0)
    
    Image_address=[datasetPath,'/',seqName,'/img']; % Image directory
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SeqReID tracker %%%%%%%%%%%%%%%%%%%%%%%%%%
    Trk=SeqReID_Tracker(idl0,Image_address,param); % SeqReID Tracker
    
    [itlf_0,fp_idl_0] = Reformat_Trajectories(Trk,idl0); % Reformating trajectories
    % to be consistent with the evaluation format used in [Dicle et al. ICCV 2013]
  
    fprintf(['Computing MOT metrics for SeqReID on ',seqName,' dataset \n']);
    mot{ii,1} = smot_clear_mot_fp(itl0,itlf_0,fp_idl_0,param.mota_th);
    fprintf('\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SMOT tracker %%%%%%%%%%%%%%%%%%%%%%%%%%%
    SMOT_Tracker
        
    if param.vis
        showitl(itlf{t,n},seqPath,'tail',5);
        % the below line is useful if you want to save the output
        %                 showitl(itlf{t,n},seqPath,'tail',5,'saveoutput',savePath);
    end
    
    fprintf(['Computing MOT metrics for SMOT on ',seqName,' dataset \n'])
    mot{ii,2} = smot_clear_mot_fp(itl0,itlf{t,n},fp_idl{t,n},param.mota_th);
    fprintf('\n');

end
%% Results (Table)
for ii = 1:length(dataSets)
    
    seqName = dataSets{ii}; % Dataset name
    fprintf(['Computing MOT metrics for SeqReID on ',seqName,' dataset \n']);
    fprintf('mme:%3g\tfn:%3g\tfp:%3g\tmota:%0.4f\n',mot{ii,1}.mme,mot{ii,1}.fn,mot{ii,1}.fp,mot{ii,1}.mota);
   
    fprintf(['Computing MOT metrics for SMOT on ',seqName,' dataset \n'])
    fprintf('mme:%3g\tfn:%3g\tfp:%3g\tmota:%0.4f\n',mot{ii,2}.mme,mot{ii,2}.fn,mot{ii,2}.fp,mot{ii,2}.mota);
    fprintf('\n');
end