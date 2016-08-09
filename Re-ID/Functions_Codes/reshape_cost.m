function JC=reshape_cost(Pr,Ix,C)
% This function re-generate and re-ordering assignment cost matrix based on 
% the marginal probabilities.
% Since, these margial probabilites are approximated, it may include
% many zeros which cannot be helpful for the ranking.  
% This function heuristically generates a sensible cost matrix such that
% it guarantees the cost values has the same order as the inverse of 
% non-zeros marginal probabilites and keep the same order as the original
% cost values for the zeros probabilites. 

N_t=size(Pr,2);
JC=C;

for i=1:N_t
    Ci=C(i,Ix{i});
    Cj=C(i,setdiff(1:size(C(i,:),2),Ix{i}));
    mvli=min(Ci);
    if ~isempty(Cj)
    mvlj=min(Cj);
    end
    
    [Vl,inx]=sort(Ci);
    [Vlp,inxp]=sort(Pr{i},'descend');
    inxx=inxp(Vlp~=0);
    O_inx=Ix{i}(inx);
    J_inx=Ix{i}(inxx);  
    indx=[J_inx setdiff(O_inx,J_inx,'stable')];
    sizz=size(indx,1);
    if ~isempty(Cj)
    JC(i,indx)=mvli:(mvlj-eps-mvli)/(sizz-1):mvlj-eps;
    else
        JC(i,indx)=Vl;
    end
        
end
