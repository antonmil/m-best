function F_Pr=mBest_Marginal_Probabilty_Calculator(M,N)

N_T=size(M,1); % Number of query images
F_Pr=cell(1,N_T);% Final Probabilities
msp=M{1,1}.Meas_edge(2,end); % number of k
Aeq=[]; beq=[]; % Equality constraints Aeq*X=beq
f=[];% costs for edges
Edge_inx=zeros(1,N_T); % Accumulator for the number of edges 
Edge_vec0=cell(1,N_T);
for nn=1:N_T
    Aeq=blkdiag(Aeq,M{nn}.A_Eq_Const); % Equality constraint
    beq=[beq;M{nn}.b_Eq_Const]; %#ok<*AGROW> % Equality constraint
    Edge_inx(nn)=size(f,1);
    f=[f;M{nn}.Costs]; % Linear Objective
    ind0=find(M{nn}.Meas_edge(2,:)==1);
    Edge_vec0{nn}=Edge_inx(nn)+ind0;
    F_Pr{1,nn}=zeros(size(ind0,2),1);
end

N_E=size(f,1);% Total number of edges
A=[];b=[]; % Inequality constrains A*X<=b
if N_T>1
    for ss=1:msp
        T_Meas_Set=[];
        for nn=1:N_T
            T_Meas_Set=[T_Meas_Set M{nn}.Meas_edge(1,M{nn}.Meas_edge(2,:)==ss)];
        end
        T_Meas_Set=unique(T_Meas_Set);
        T_Meas_Set_wo_zero=T_Meas_Set(~(T_Meas_Set==0));
        sizTM=size(T_Meas_Set_wo_zero,2);
        for mm=1: sizTM
            Edge_vec=cell(1,N_T);
            for nn=1:N_T
                ind=find(M{nn}.Meas_edge(2,:)==ss&M{nn}.Meas_edge(1,:)==T_Meas_Set_wo_zero(mm));
                if ~isempty(ind)
                    Edge_vec{nn}=Edge_inx(nn)+ind;
                end
            end
            Edge_vec(cellfun(@isempty,Edge_vec))=[];
            N_T_C=size(Edge_vec,2);
            if N_T_C>1
                Edge_comb=combvec(Edge_vec{:});
                Num_ineq=size(Edge_comb,2);
                b0=ones(Num_ineq,1);
                Edge_comb_y = reshape(Edge_comb',1,Num_ineq*N_T_C);
                Edge_comb_x=repmat(1:Num_ineq,[1 N_T_C]);
                A0=sparse(Edge_comb_x,Edge_comb_y,ones(size(Edge_comb_y)),Num_ineq,N_E);
                A=[A;A0]; % Inequality constrain
                b=[b;b0]; % Inequality constrain
            end
        end
        
    end
end

options = optimoptions('intlinprog','Display','off','CutGeneration','none','BranchingRule','mostfractional');
[candidates, values] = BinIntMBest(f,A,b,Aeq,beq,options,N); % Calculating
% m-best solutions of the linear binary program with one-to-one constraints

%% Calculating marginal probabilites for query images
for hc=1:length(values)
    x2 = candidates(:,hc);
    f_prr=exp(-x2'*f);
    
    for nn=1:N_T
     F_Pr{1,nn}(logical(x2(Edge_vec0{nn})),1)=F_Pr{1,nn}(logical(x2(Edge_vec0{nn})),1)+f_prr;
    end
end

F_Pr=cellfun(@(x) x/sum(x), F_Pr, 'UniformOutput', false); 



