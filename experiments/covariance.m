
% covariance experiments


x = rand(1,1000) * 3 + 5;
mean_x = mean(x);
std_x  = std(x);
y = rand(1,1000) * -4 + 0.3 * x;
mean_y = mean(y);
std_y  = std(y);
figure(3);
plot(x,y,'.');


x_fin = fin_res(1).dat(1,:);
y_fin = fin_res(1).dat(2,:);
figure(4);
plot(x_fin, y_fin, '.');
mean(x_fin)
mean(y_fin)
std(x_fin)
std(y_fin)
cov([x_fin; y_fin]')







res.data_ratio  = sum(sum(abs(p * dat))) / sum(sum(abs(p2 * dat)));

%   Cholesky-like decomposition for covariance matrix.
%   T = CHOLCOV(SIGMA) computes T such that SIGMA = T'*T.
% T = cholcov(cov(dat));


% bootstrapping loop
for i = 1:n_iter
    
    %create test_space with same statisical properties as given data and compute ratio
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
    res.test_ratios(i)   = sum(sum(abs(p * dat_test))) / sum(sum(abs(p2 * dat_test)));
    
    %shuffle data and compute ratio      --> old all over shuffeling
    %   dat_shuf             = reshape(dat(randperm(length(dat(:)))), size(dat) );  %

    % das gab die nicht schönen ergebnisse
     dat_shuf = NaN(size(dat));
%     for j = 1:size(dat,1)
%         dat_shuf(j,:) = dat(j,randperm(size(dat,2)));
%     end
%     res.shuf_ratios(i)  = sum(sum(abs(p * dat_shuf))) / sum(sum(abs(p2 * dat_shuf)));


% jetzt noch mal versuchen, wenn ich zwischen 1 und 8 alle vertausche und
% zwischen dem rest auch
    fl18 = dat([1 8],:);
    fl18 = fl18(:);
    dat_shuf([1 8],:) = reshape(fl18(randperm(length(fl18))),2,148);
    fl18 = dat([2:7 9:11],:);
    fl18 = fl18(:);
    dat_shuf([2:7 9:11],:) = reshape(fl18(randperm(length(fl18))),9,148);
    res.shuf_ratios(i)  = sum(sum(abs(p * dat_shuf))) / sum(sum(abs(p2 * dat_shuf)));

    %create random data and compute ratio
    dat_rand            = rand(size(dat));
    res.rand_ratios(i)  = sum(sum(abs(p * dat_rand))) / sum(sum(abs(p2 * dat_rand)));
    
    % create data with same statistical properties, also second order
    if 0
        % mu = dat;
        mu = ones(size(dat,1),1) * mean(dat);
        
        %jetzt versuche ich mal aus einer verteilung zu ziehen, die unseren daten
        %ähnelt
        %     x_ = 0:0.000001:1;
        %     p_ = mle(dat(:), 'dist','tlocationscale', 'alpha',0.05);  % Fit t location-scale distribution
        %     dist = pdf('tlocationscale',x_,p_(1), p_(2), p_(3));
        %     perms = randperm(length(dist));
        %     dat_rcov = reshape(dist(perms(1:size(dat,1)* size(T,1))), size(dat,1), size(T,1)) * T;
        
        
        dat_rcov             = randn(size(dat,1),size(T,1)) * T + mu;
        % rescale it to the same range as dat
        % NOTE, is this correct ?
        dat_rcov = dat_rcov - min(min(dat_rcov));
        dat_rcov = dat_rcov ./ max(max(dat_rcov));
        res.rcov_ratios(i)   = sum(sum(abs(p * dat_rcov))) / sum(sum(abs(p2 * dat_rcov)));
    else
        dat_bla             = reshape(dat(randperm(length(dat(:)))), size(dat) );  %
        res.rcov_ratios(i)   = sum(sum(abs(p * dat_bla))) / sum(sum(abs(p2 * dat_bla)));
    end
end