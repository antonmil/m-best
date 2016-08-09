function Demo_ReID_mbst(varargin)

%% Code for Person Re-Identification (ReID) using m-best solutions  

% It requires:
% 1 - Assinment matrix representing a score (or cost) between all query (
%     rows in the matrix) and gallery (columns of the matrix) images. 
% 2 - Number of solutions (mbest parameter) 
% 3 - Gurobi ILP solver

% The code is tested under Windows 10 (64bit)

% If you use this code for your research, please, cite:
% S. H. Rezatofighi, A. Milan, Z. Zhang, Q. Shi, A. Dick, I. Reid, "Joint 
% Probabilistic Matching Using m-Best Solutions ", IEEE Conference 
% on Computer Vision and Pattern Recognition (CVPR), 2016.

% Author: S. Hamid Rezatofighi
% Last updated: July 18, 2016

% Note: This code package can be used for person ReID as an additional step 
% to enhance the matching results using joint information. 
% This version of code has been implemented for the case when all query
% persons exist in the gallery images. Therefore, the input should be a 
% square assignment matrix. However, its extension to the other cases is 
% straightforward, e.g. by adding a dummy node to the partitie 
% graph to deal with non-square assignment matrix cases, or augment the 
% bipartite graph to multipartite when there a multiple instances of same 
% person in the gallery set. For more information, you can check the
% following paper

% S. H. Rezatofighi, A. Milan, Z. Zhang, Q. Shi, A. Dick, I. Reid, "Joint 
% Probabilistic Data Association Revisited", IEEE International Conference 
% on Computer Vision (ICCV), 2015.

% For questions contact the author at: hamid.rezatofighi@adelaide.edu.au

AddPath; % Add the necessary paths 

narginchk(0, 2)
if ~nargin||isempty(varargin{1})
    BDs = loadTable('cuhk03', 1);% Loading assignment matrix (as cost)
else
    BDs = varargin{1};% Assignment matrix (as cost)
end
if ~nargin||isempty(varargin{2})
    mbest=100; % Number of m-best solutions, it should be m>>1 
else
    mbest = varargin{2};% Number of m-best solutions, it should be m>>1 
end

[BhD,Tar_inx] = Objective_Constraints (BDs);
Final_probabilty=mBest_Marginal_Probabilty_Calculator(BhD,mbest);
JBD=reshape_cost(Final_probabilty,Tar_inx,BDs);

cmc = evaluate_pwdist(BDs');
mbst_cmc = evaluate_pwdist(JBD');

plot(1:size(cmc,1),cmc,'*-k','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(1:size(mbst_cmc,1),mbst_cmc,'d-g','LineWidth',2,'MarkerFaceColor','r','MarkerSize',5)
legend('Original Cost','mbst on Original Cost','Location','SouthEast')
xlabel('Rank')
ylabel('Recognition rate (%)')
title ('CMC curve')

end