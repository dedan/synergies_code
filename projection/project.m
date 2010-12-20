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

n_rand      = 10000;
test_data   = 0;

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
res.data_ratio  = sum(sum((p * dat).^2)) / sum(sum((p2 * dat).^2));
res.means       = mean(dat,2);



% bootstrapping loop
for i = 1:n_iter
    
    %create test_space with same statisical properties as given data and compute ratio
    if test_data
        x        = linsolve(basis,dat);
        mean_x   = mean(x,2);
        std_x    = std(x,0,2);
        x_basis  = NaN(size(x));
        for j = 1:size(x,1)
            distribution = (randn(1,n_rand) * std_x(j)) + mean_x(j);
            perms        = randperm(n_rand);
            x_basis(j,:) = distribution(perms(1:size(x,2)));
        end
        dat_test         = basis * x_basis;
        
        % put noise on test data, otherwise very high value
        noise                = randn(size(dat_test)) * mean(dat_test(:)) * noise_value;
        dat_test             = dat_test + noise;
        res.test_ratios(i)   = sum(sum((p * dat_test).^2)) / sum(sum((p2 * dat_test).^2));
    else
        res.test_ratios(i)   = 0;
    end
    
    % random with same overall mean and variance
    dat_rand    = normrnd(mean(dat(:)), std(dat(:)), size(dat,1), size(dat,2));
    res.rand_ratios(i)  = sum(sum((p * dat_rand).^2)) / sum(sum((p2 * dat_rand).^2));

    % random with same distribution for every channel
    dat_rand = NaN(size(dat));
    for j = 1:size(dat,1)
        dat_rand(j,:)   = normrnd(mean(dat(j,:)), std(dat(j,:)), 1, size(dat,2));
    end
    res.test_ratios(i)  = sum(sum((p * dat_rand).^2)) / sum(sum((p2 * dat_rand).^2));
    
    % shuffled data
    dat_shuf            = reshape(dat(randperm(length(dat(:)))), size(dat) );
    res.shuf_ratios(i)  = sum(sum((p * dat_shuf).^2)) / sum(sum((p2 * dat_shuf).^2));
    res.shuf_means(:,i) = mean(dat_shuf,2);
    
    % shuffle only between sessions
    dat_rows  = NaN(size(dat));
    for j = 1:size(dat,1)
        dat_rows(j,:)   = dat(j,randperm(size(dat,2)));
    end
    res.rcov_ratios(i)  = sum(sum((p * dat_rows).^2)) / sum(sum((p2 * dat_rows).^2));
    res.rcov_means(:,i) = mean(dat_rows,2);
end



