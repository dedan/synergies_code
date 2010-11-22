function result = nmf_explore(dat, config)

result = struct;
col = struct;

SS = sum(sum(dat));

for j = 1:config.Niter_exploration;
    
    % first run the multiplicative nnmf alorithm mit different
    % initializations, always only for 5 iterations in order to
    % explore the search space
    [W0,H0] = nnmf(dat',config.dimensions,'replicates',100,'options',config.opt,'algorithm','mult');
    
    % take best result of exploration episodes and converge using
    % the alternating update rule
    [W H]  = nnmf(dat', config.dimensions,'w0',W0,'h0',H0,'algorithm','als');
    col(j).syns = W';
    col(j).errs = (sum(sum(abs(dat' - W * H))) / SS) * 100;
end

% sort the error values and choose twenty best results
result.errs = [col.errs];
[~, idx]    = sort(result.errs);
best        = col(idx(1:config.n_best));

% group, plot and save the results  --> show nmf stability
grouped     = group(best, 'syns');
result.flat = flatten(grouped);
result.syns = grouped(1).center;
result.std  = std(hist(grouped(1).idx, config.dimensions));


