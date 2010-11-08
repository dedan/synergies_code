function res = shuffle(mat)

    % randomly shuffle a matrix

    flat = mat(:);
    res = reshape(flat(randperm(length(flat))),size(mat,1),size(mat,2));
end