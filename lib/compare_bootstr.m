function dist = compare_bootstr(dat, n_iter)

dat = shuffle_inc(dat);
for i = 1:n_iter
    rands   = randperm(size(dat,1));
    syns    = normr(dat(rands(1:2),:));
    dist(i,:) = matchNscore(syns(1,:)',syns(2,:)');    
end
    
    