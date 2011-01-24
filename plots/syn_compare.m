function idx = syn_compare(syn1, syn2, name1, name2)

% pairwise comparison between two groups of vectors.
% second group will be matched to sorting of first group.

% returns indices how to sort the second matrice, according to the first.
% I use it, as I know that the pca results are already sorted, to sort the
% nmf results according to them

% similarity measured by dot product of normed vectors

syn1 = normr(syn1);
syn2 = normr(syn2);

if(size(syn1) == size(syn2))
    n = size(syn1,1);
else
    error('the two groups must have the same size');
end

[syn1, syn2, scores] = match_syns(syn1, syn2);

for j = 1:n
    subplot(n+1,1,j)
    bar( [syn1(j,:)' syn2(j,:)']);

    % set y axis limits depending on negative or non-negative data
    if size(find(syn1 < 0),1) > 0 || size(find(syn2 < 0),1) > 0
        ylim([-1 1]);
    else
        ylim([0 1]);
    end
    
    title(['synergy #' int2str(j) ' matching score: ' num2str(scores(j))]);
end

subplot(n+1,1,n+1)
plot(syn1(:), syn2(:), '.');
[r p] = corrcoef([syn1(:) syn2(:)]);

title([name1 ' vs. ' name2 ', r: ' num2str(r(1,2)) ' - p: ' num2str(p(1,2))]);
xlabel(name1);
ylabel(name2);
