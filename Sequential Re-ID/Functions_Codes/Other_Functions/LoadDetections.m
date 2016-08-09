clear idl0 itl0
% Load the annotated detections
if strcmp(seqName,'AFL_S1')||strcmp(seqName,'PETS')
    load([datasetPath,'/',seqName,'/detection_gt/GT_Detections.mat'])
    itl0=loaditl([datasetPath,'/',seqName,'/detection_gt/GT_Detections.itl']);
    seqPath = [datasetPath '/' seqName '/img'];
    for k=1:length(gtInfo.frameNums)
        idl0(k).rect(:,1)=(gtInfo.X(k,(gtInfo.X(k,:)~=0))-gtInfo.W(k,(gtInfo.W(k,:)~=0))/2)';
        idl0(k).rect(:,2)=(gtInfo.Y(k,(gtInfo.Y(k,:)~=0))-gtInfo.H(k,(gtInfo.H(k,:)~=0)))';
        idl0(k).rect(:,3)=(gtInfo.W(k,(gtInfo.W(k,:)~=0)))';
        idl0(k).rect(:,4)=(gtInfo.H(k,(gtInfo.H(k,:)~=0)))';
        
        idl0(k).xy(:,1)=(gtInfo.X(k,(gtInfo.X(k,:)~=0)))';
        idl0(k).xy(:,2)=(gtInfo.Y(k,(gtInfo.Y(k,:)~=0))-gtInfo.H(k,(gtInfo.H(k,:)~=0))/2)';
    end
elseif strcmp(seqName,'AFL_S2')||strcmp(seqName,'AFL_S3')
    load([datasetPath,'/',seqName,'/detection_gt/GT_Detections.mat'])
    itl0=loaditl([datasetPath,'/',seqName,'/detection_gt/GT_Detections.itl']);
    seqPath = [datasetPath '/' seqName '/img'];
else
    initialize_smot;
end