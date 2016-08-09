function cmc = evaluate_pwdist(pwsim)
% evaluate the CMC performance of pair-wise distance
% assume gallery in dim1, query in dim2
% 
N = size(pwsim, 3);
gsize = size(pwsim, 1);    
cmc = zeros(gsize, N);
for i = 1:N
    [~, order] = sort(pwsim(:, :, i));
    match = (order == repmat(1:gsize, [gsize, 1]));
    cmc(:, i) = cumsum(sum(match, 2)./gsize);
end