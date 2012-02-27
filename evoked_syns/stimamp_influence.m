% % investigate the influence of the stimulations strength

% we were regularly asked whether we are sure that we stimulate only
% single neurons. Yifat thinks that we stimulate single neurons
% or at least small groups. This is an important point because if we would
% already stimulate groups of neurons the synergies we obsever could
% arise from (linear) combinations of activity patterns of neurons.
% To check this the influence of the stimulation amplitude is investigated.
% -> The larger the amplitude the more neurons should be stimulated
% -> if the activation pattern is a combination of singel neuron patterns, it should
%    change when more neurons are activated.
% -> therefore this function checks whether the stimulation amplitude influences
%    the direction of the activation pattern (evidence that it is combined from single
%     neuron activation patterns) or that it influences the only the amplitude
% the idea is to use a svd to check the rank of the matrix which indicates the linear
% dependency of its rows

%  calling example
%  stimamp_influence('~/projects/yifat_paper/results/data/', ...
%                    '~/projects/yifat_paper/results/evoked_syns/',  {'vega', 'chalva'})
function stimamp_influence(inpath, outpath, monks)

    for i = 1:length(monks)
        load([inpath 'allevoked_data_' monks{i}]);

        % a session can have different locations.
        % for each location we want to compare responses for different amplitudes
        res = read_resps(resps);

        % combine all responses in one matrix and shuffle it
        conc = vertcat(res.ry);
        shuffler = randperm(size(conc, 1));
        conc = conc(shuffler, :);
        ind = 1;

        % for all locations
        for j = 1:length(res)

            % check linear dependency
            s = svd(res(j).ry);
            eigenvalue_fraction(j) = s(1) / sum(s);

            % and also for random activation patterns
            l = size(res(j).ry, 1);
            ry_shuff = conc(ind:ind+l-1, :);
            ind = ind + l;
            s = svd(ry_shuff);
            shuff_eigenvalue_fraction(j) = s(1) / sum(s);
        end

        f = figure('Visible', 'off');
        subplot 211
        hist(eigenvalue_fraction, 100)
        title('distribution of first eigenvalue fraction for locations')
        subplot 212
        hist(shuff_eigenvalue_fraction, 100)
        title('same for random activation patterns')
        saveas(f, [outpath 'stimamp_influence_' monks{i} '.png'])
        close(f)

        disp([monks{i} ' --> ranksum test original and shuffled distribution'])
        [p,h] = ranksum(eigenvalue_fraction, shuff_eigenvalue_fraction);
        disp(sprintf('\tp: %f - h: %f\n\n', [p, h]))
    end
end


% split responses for the different locations
function res = read_resps( resps )

    ampx    = [25 50:50:250];
    N       = length(ampx);
    reci    = 1;
    curindx = 1;
    res = struct;
    while curindx <= length(resps),
        sessindx = resps(curindx).id;
        depth    = resps(curindx).location.depth;
        el       = resps(curindx).electrode;

        indxs    = get_indx(resps, curindx, sessindx, depth, el);
        res(reci).rx = [];
        res(reci).ry = [];
        for i=1:length(indxs),
            ipos = find(ampx == abs(resps(indxs(i)).amp));
            res(reci).rx = [res(reci).rx; ampx(ipos)];
            res(reci).ry = [res(reci).ry; resps(indxs(i)).response];
        end
        reci = reci+1;
        curindx = indxs(end)+1;
    end
end


% get all indeces of same session, depth, electrode
function indxs = get_indx( resps, curi, sessindx, depth, el)

    indxs(1) = curi;
    curi = curi+1;
    while curi <= length(resps) && resps(curi).electrode == el && ...
     resps(curi).id == sessindx &&  resps(curi).location.depth == depth,
        indxs(end+1) = curi;
        curi = curi+1;
    end
end
