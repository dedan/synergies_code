function res = project(dat, basis, n_iter, noise_value)             % takes already column vectors

% original data is compared to different shuffeling, test data and random
% data.

% if random data has the same mean and variance as dat, it has the same
% distribution as dat when shuffled without restriction

% random data which has same mean in variance over the channels has the
% same distribution as da shuffled only within channels

% test data is data created by linear combinations of the basis vectors
% with similar statistical properties as dat. This can give an estimate of
% the upper limit of the ratio, when different amounts of gaussian noise
% are added on the test data

dat = dat';
basis = basis';

if size(basis,1) ~= size(dat,1)
    error('first dimension has to be the same');
end

% orthonormalization and creation of orthogonal complement
basis    = orth(basis);
ortcmp   = null(basis');

% create projections 
p        = basis * basis';
p2       = ortcmp * ortcmp';

% compute the ratio for the data
res.ratio_dist  = sum((p * dat).^2) ./ sum((p2 * dat).^2);
res.dat_ratio   = sum(sum((p * dat).^2)) / sum(sum((p2 * dat).^2));
    
% bootstrapping loop
for i = 1:n_iter
        
    % random with same overall mean and variance
    dat_rand            = normrnd(mean(dat(:)), std(dat(:)), size(dat,1), size(dat,2));
    res.rand_ratios(i)  = sum(sum((p * dat_rand).^2)) / sum(sum((p2 * dat_rand).^2));

    % random with same distribution for every channel
    dat_rand = NaN(size(dat));
    for j = 1:size(dat,1)
        dat_rand(j,:) = normrnd(mean(dat(j,:)), std(dat(j,:)), 1, size(dat,2));
    end
    res.rand_chan_ratios(i) = sum(sum((p * dat_rand).^2)) / sum(sum((p2 * dat_rand).^2));
    
    % shuffled data
    dat_shuf                = reshape(dat(randperm(length(dat(:)))), size(dat) );    
    res.shuf_ratios(i)      = sum(sum((p * dat_shuf).^2)) / sum(sum((p2 * dat_shuf).^2));
    
    % shuffle only in channel
    dat_rows  = NaN(size(dat));
    for j = 1:size(dat,1)
        dat_rows(j,:) = dat(j,randperm(size(dat,2)));
    end
    res.shuf_chan_ratios(i)  = sum(sum((p * dat_rows).^2)) / sum(sum((p2 * dat_rows).^2));
end



