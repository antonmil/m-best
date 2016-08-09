function [BhD,Tar_inx] = Objective_Constraints (BDs)       
nqury=size(BDs,2); % number of query images
ngllry=size(BDs,2); % number of gallery images
if nqury~=ngllry
    error('Current implementation can handle a square assignment matrix only')
    % The function mBest_Marginal_Probabilty_Calculator need to be modified 
    % to deal with dummy node for non-square assignment matrix
end
BD_g=(BDs<=inf); % Check for Inf values (no possible connection between the
% nodes) in cost. 

BhD=cell(nqury,1);
Tar_inx=cell(1,nqury);


for i=1:nqury
    Tar_inx{i}=find(BD_g(i,:)==1); % Indecies of possible gallery images  
    % for the i-th query image
    BhD{i}.Costs=BDs(i,Tar_inx{i})';
    BhD{i}.Meas_edge=[Tar_inx{i};ones(1,size(Tar_inx{i},2))];
    BhD{i}.A_Eq_Const=sparse(ones(1,size(Tar_inx{i},2)));
    BhD{i}.b_Eq_Const=1;
end