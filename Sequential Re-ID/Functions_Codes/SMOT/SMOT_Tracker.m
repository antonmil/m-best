
n=N;t = TRIALS;

fprintf('Processing Dataset:%s  Method:%s  Noise:%0.2f|%0.2f  Trial:%d\n',...
    seqName,method,fp(n),fn(n),t);
fprintf('---------------------------------------\n');
if  n==1 && t>1
    % noiseless case is same for all trials. copy.
    itlf{t,n} = itlf{1,n};
    etime(t,n) = etime(1,n);
    mot(t,n) = mot(1,n);
else
    noise(n).fp = fp(n);
    noise(n).fn = fn(n);
    noise(n).gn = 0;
    
    [idlp,fp_idl{t,n},fn_idl{t,n}] = idladdnoise(idl0,noise(n));
    
    
    % do the stitching
    if isequal(method, 'ksp')
        % initialize ksp parameters
        %         initialize_ksp;
        [itlf{t,n} etime(t,n)] = ksp_associate(idlp,param);
    else
        [itlf{t,n} etime(t,n)] = smot_associate(idlp,param);
    end
    
    
    
    
    fprintf('Process time: %g \n\n',etime(t,n));
 
end

