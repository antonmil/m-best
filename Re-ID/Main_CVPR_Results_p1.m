%% Code for Person Re-Identification (ReID) using m-best solutions  

% It requires:
% 1 - Assignment matrix representing a score (or cost) between all query (
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
% straightforward, e.g. by adding a dummy node to the partite 
% graph to deal with non-square assignment matrix cases, or augment the 
% bipartite graph to multipartite when there a multiple instances of same 
% person in the gallery set. For more information, you can check the
% following paper

% S. H. Rezatofighi, A. Milan, Z. Zhang, Q. Shi, A. Dick, I. Reid, "Joint 
% Probabilistic Data Association Revisited", IEEE International Conference 
% on Computer Vision (ICCV), 2015.

% For questions contact the author at: hamid.rezatofighi@adelaide.edu.au

%% Demo

% For a test, try an example assignment cost "BDs" for ReID and use the 
% following functions 
% -------------------------------------------------------------------------
% [BhD,Tar_inx] = Objective_Constraints (BDs);
% Final_probabilty=mBest_Marginal_Probabilty_Calculator(BhD,mbest);
% JBD=reshape_cost(Final_probabilty,Tar_inx,BDs);

% and plot CMC curves before and after marginalization using m-best
% solutions

% cmc = evaluate_pwdist(BDs');
% mbst_cmc = evaluate_pwdist(JBD');

% plot(1:size(cmc,1),cmc,'*-k','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
% hold on
% plot(1:size(mbst_cmc,1),mbst_cmc,'d-g','LineWidth',2,'MarkerFaceColor','g','MarkerSize',5)
% legend('Original Cost','mbst on Original Cost','Location','SouthEast')
% xlabel('Rank')
% ylabel('Recognition rate (%)')
% title ('CMC curve')

% -------------------------------------------------------------------------
%% CVPR Results for iLIDS, 3DPeS, VIPeR, CUHK01, CUHK03
close all
clear all
clc

AddPath; % Add the necessary paths 

DATASET = {'ilids', '3dpes', 'viper','cuhk01', 'cuhk03'}; % Dataset name
d_siz=size(DATASET,2);
cmc=cell(1,d_siz);
mbst_cmc=cell(1,d_siz);

No_Ex=1:10;% No of experiments, similar to [Paisitkriangkrai et al., CVPR 2015]
p_time=zeros(d_siz,length(No_Ex)); % Pre-allocation for time
mbest=100; % Number of m-best solutions, it should be m>>1

for jd=1:d_siz
    for trialNo=No_Ex
        disp([DATASET{jd},num2str(trialNo)])
        clear Final_probabilty BDs BD_g BhD Tar_inx 
        BDs = loadTable(DATASET{jd}, trialNo); % Loading assignment matrix (as cost)
        % Exactly the same cost matrix as [Paisitkriangkrai et al., CVPR 2015]
        
        [BhD,Tar_inx] = Objective_Constraints (BDs); % Building the linear objective 
        % and one-to-one constraints for linear binary programming  

        tic;
        Final_probabilty=mBest_Marginal_Probabilty_Calculator(BhD,mbest);
        % Calculating marginal probabilities using m-best solutions
        
        JBD=reshape_cost(Final_probabilty,Tar_inx,BDs); % Re-generating and 
        % re-ordering assignment cost matrix based on marginal probabilities.

        p_time(jd,trialNo)=toc;
        
        cmc{jd}(:,trialNo) = evaluate_pwdist(BDs'); % CMC curve before 
        % marginalization
        
        mbst_cmc{jd}(:,trialNo) = evaluate_pwdist(JBD');% CMC curve after 
        % marginalization
        
        
        fprintf('Rank-1 Accuracy for original costs: %f\n', cmc{jd}(1,trialNo));
        fprintf('Rank-1 Accuracy for m-best marginal costs: %f\n', mbst_cmc{jd}(1,trialNo));
    end
    save([pwd,filesep,'Results',filesep,'All_ReID_p1_m',num2str(mbest),'.mat'],...
    'cmc','mbst_cmc','mbest','DATASET','p_time')
end

%% Results and dispaly
rank_v=[1 2 5 10 20];
for jd=1:d_siz
    Og=100*mean(cmc{jd},2);
    Jm=100*mean(mbst_cmc{jd},2);
        disp(DATASET{jd})
        
        hfig=figure;
        
        plot(1:size(Og,1),Og,'*-k','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
        hold on
        
        plot(1:size(Jm,1),Jm,'d-g','LineWidth',2,'MarkerFaceColor','g','MarkerSize',5)
        
        set(gca,'FontName','Times','FontSize',20),
        set(gca, 'box', 'off')
        
        legend('Original Cost','mbst on Original Cost','Location','SouthEast')
        % legend('boxoff')
        xlabel('Rank')
        ylabel('Recognition rate (%)')
        grid on
        axis([1 25 0 101])
        set(gca,'YTick',0:10:100),
        set(gca,'XTick',0:5:25),
       
%             print(hfig,'-dpsc',[pwd,'\Results\',DATASET{jd},'.eps'])
  
            title(DATASET{jd})
            printResult(Og, Jm, DATASET{jd},mean(p_time(jd,:)),rank_v)

end

