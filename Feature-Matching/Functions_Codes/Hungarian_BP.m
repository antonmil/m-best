function assign = Hungarian_BP(K, Ct, asgT, options)
[NofNodes, ~] = size(Ct);
[~, NofStates] = size(Ct);
NodeClusters = 0:1:(NofNodes - 1);
NodeClusters = mat2cell(NodeClusters, 1, ones(1, NofNodes)); %#ok<MMTC>
Npfunc = full(diag(K));
Npfunc2 = reshape(Npfunc, [NofStates, NofNodes]);
NodePotentials2 = mat2cell(Npfunc2, NofStates, ones(1, NofNodes));
EdgesMat = zeros(NofNodes, NofNodes);
K1 = K;

K = K - diag(diag(K));
ecnt2 = 1;
Edges2 = {};

EdgePotentials2 = {};

for i=1:NofNodes
    for j=1:NofNodes
        if(i==j)
            continue;
        end
        efpunc2 = K(((i - 1) * NofStates) + (1:NofStates), ((j - 1) * NofStates) + (1:NofStates) )';
        if(nnz(efpunc2) ~= 0)
            efpunc2 = efpunc2 - 100 * eye(NofStates, NofStates);

            Edges2{ecnt2} = [i - 1, j - 1];
            efpunc2 = full(efpunc2(:));
            EdgePotentials2{ecnt2} = efpunc2(:);
            EdgesMat(i,j) = 1;
            EdgesMat(j, i) = 1;
            ecnt2 = ecnt2 + 1;
        end
    end
end






NodeStates = NofStates * ones(1, NofNodes);
Factors2 = [ NodeClusters, Edges2];%, Triples];


Potentials2 = [NodePotentials2, EdgePotentials2];%, TPFunctions];%, TPfuncs];



%[decode1, ~, ~, dual1] = MatchingBP(NofNodes, NodeStates, Factors1, Potentials1, options.outIter, options.innerIter,  -1e20);

decode = BaBMatchingSolver(NofNodes, Factors2,Potentials2,[],[],options);
%[decode2, ~, ~, dual2] = MatchingBP(NofNodes, NodeStates, Factors2, Potentials2, options.outIter, options.innerIter,  -1e20);





assign.alg='gmBP';
assign.X = zeros(NofStates, NofStates);
for i=1:NofNodes
    assign.X(decode(i) + 1, i) = 1;
end
X=assign.X;




acc = matchAsg(assign.X , asgT);
assign.acc = acc;
assign.obj = X(:)' * K1 * X(:);


end
