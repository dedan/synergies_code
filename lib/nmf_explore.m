function result = nmf_explore(dat, conf)

result = struct;
col    = struct;

centered = dat - repmat(mean(dat), size(dat, 1), 1);
SS = sum(sum(centered.^2));

for j = 1:conf.Niter_exploration;

    % first run the multiplicative nnmf alorithm mit different
    % initializations, always only for the amount of iterations given in opt iterations in
    % order to explore the search space
    [W0,H0]     = nnmf(dat',conf.dim,'replicates',100,'options',conf.opt,'algorithm','mult');

    % take best result of exploration episodes and converge using
    % the alternating update rule
    [W H]       = nnmf(dat', conf.dim,'w0',W0,'h0',H0,'algorithm','als');
    col(j).syns = W';
    err         = dat' - W * H;
    col(j).errs = (sum(sum(err.^2)) / SS) * 100;
end

% sort the error values and choose twenty best results
result.errs = [col.errs];
[~, idx]    = sort(result.errs);
best        = col(idx(1:conf.n_best));

% group, plot and save the results  --> show nmf stability
grouped     = group(best, 'syns');
result.flat = vertcat(grouped.dat);
result.syns = grouped(1).center;
result.std  = std(hist(grouped(1).idx, conf.dim));


