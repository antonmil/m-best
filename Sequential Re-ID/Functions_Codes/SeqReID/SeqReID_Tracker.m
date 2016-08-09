function Trk=SeqReID_Tracker(idl0,Image_address,param)

Gate = param.G; % Mamximum distance to consider matching candidates (pixels)
term_f = param.term_f; % Termination parameter
cost0 = param.min_dist; % Distance to missed detection
Frme =size(idl0,2); % Number of frames

colorord = distinguishable_colors(1000); % Generate distinguishable colors
% for different targets

I = LoadImg(Image_address,1); % Load each Image frame

Feat=cell(1,Frme);
XY=cell(1,Frme);
[Feat{1,1},XY{1,1}] = HSV_HOG(I,idl0,1); % Calculate HSV_HOG features
% and extract bounding box positions

% Pre-allocation
MTch_Tg=cell(1,Frme);
F_Pr=cell(1,Frme);
Ntg=size(Feat{1,1},1); % Number of targets
Trk=cell(1,Ntg);
Exist_trg=1:Ntg; % Existing track IDs
Oc_cnt=zeros(1,Ntg); 


for i=Exist_trg
    Trk{i}=[i;0;1]; % Building Trajectories
end

if param.vis
    TrackVisualizer(I,Trk,idl0,1,colorord);
    drawnow
end

for k=2:Frme
    I = LoadImg(Image_address,k); % Load each Image frame
    [Feat{1,k},XY{1,k}] = HSV_HOG(I,idl0,k); % Calculate HSV_HOG features
    
    D_si2=size(Feat{1,k},1);
    MTch_Tg{k}=cell(max(Exist_trg),1);
    Mes_Tar=false(D_si2,max(Exist_trg));
    for i=Exist_trg
        MTch_Tg{k}{i}.Costs(1)=cost0;
        MTch_Tg{k}{i}.Meas_edge=[0;1];
        for j=1:D_si2
            
            if Trk{i}(1,end)~=0
                Feat_pst=Feat{1,Trk{i}(3,end)}(Trk{i}(1,end),:);
                xy=XY{1,Trk{i}(3,end)}(Trk{i}(1,end),:);
            else
                lnz=find(Trk{i}(1,:)~=0);
                Feat_pst=Feat{1,Trk{i}(3,lnz(end))}(Trk{i}(1,lnz(end)),:);
                xy=XY{1,Trk{i}(3,lnz(end))}(Trk{i}(1,lnz(end)),:);
            end
            if pdist2(xy,XY{1,k}(j,:))<=Gate
                MTch_Tg{k}{i}.Costs = [MTch_Tg{k}{i}.Costs;pdist2(Feat_pst,...
                    Feat{1,k}(j,:))];
                MTch_Tg{k}{i}.Meas_edge=[MTch_Tg{k}{i}.Meas_edge [j;1]];
                Mes_Tar(j,i)=true;
            end
        end
        
        MTch_Tg{k}{i}.Hypo=MTch_Tg{k}{i}.Meas_edge(1,:)';
        MTch_Tg{k}{i}.Prob=exp(-MTch_Tg{k}{i}.Costs);
        MTch_Tg{k}{i}.A_Eq_Const=sparse(ones(1,size(MTch_Tg{k}{i}.Meas_edge,2)));
        MTch_Tg{k}{i}.b_Eq_Const=1;
    end
    ixxd=(~cellfun(@isempty,MTch_Tg{k}));
    F_Pr{k}=cell(1,size(MTch_Tg{k},1));
    Mes_Tar2=Mes_Tar(:,ixxd);[Umt, Vmt]=size(Mes_Tar2);
    Mes_Tar3=[false(Vmt,Vmt+Umt);Mes_Tar2 false(Umt,Umt)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% m-best matching %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F_Pr{k}(ixxd) =Approx_Multiscan_JPDA_Probabilities(Mes_Tar3,MTch_Tg{k}(ixxd),100);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    assigned_meas=[];% Exist_trg_new=[];
    for i=Exist_trg
        [Vlu,Ixx]=max(F_Pr{1,k}{i});
        Ix=MTch_Tg{k}{i}.Hypo(Ixx);
        Trk{i}=[Trk{i} [Ix;Vlu;k]];
        if Ix~=0
            assigned_meas=[assigned_meas Ix]; %#ok<AGROW>
            Oc_cnt(i)=0;
        else
            Oc_cnt(i)=Oc_cnt(i)+1;
        end
    end
    n_trgt=setdiff(1:D_si2,assigned_meas);
    trc_siz=size(Trk,2);
    for ij=1:length(n_trgt)
        Trk{trc_siz+ij}=[n_trgt(ij);0;k];
        Oc_cnt(trc_siz+ij)=0;
    end
    Exist_trg_new=1:size(Trk,2);
    Exist_trg=Exist_trg_new(Oc_cnt<term_f);
    
    if param.vis
        TrackVisualizer(I,Trk,idl0,k,colorord);
        drawnow
    end
    
    
end

