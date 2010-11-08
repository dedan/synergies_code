function shuffled = shuffle_inc(mat)

% shuffle the matrix, but only within columns


shuffled = NaN(size(mat));

for i = 1:size(mat,2)
    perms = randperm(size(mat,1));
    shuffled(:,i) = mat(perms,i);
end