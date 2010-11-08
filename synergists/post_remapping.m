function fin_res = post_remapping(sep_results, config)

% NOTE because for vega the mapping of the muscles changed after session 25
% the results have to be separated again, from here I continue without
% these sessions
fin_res     = struct;
remap_idx   = struct;
for i = 1:length(sep_results)
    tmp = sep_results(i).resp;
    tmp = [tmp.info];
    if strcmp(config.monkey, 'vega')
        remap_idx(i).before  = find([tmp.id] < 26);
        remap_idx(i).after   = find([tmp.id] > 25);
    else
        remap_idx(i).before  = [];
        remap_idx(i).after   = [tmp.id];
    end
end


for i = 1:length(sep_results)
    if config.normalization == 1
        fin_res(i).dat  = sep_results(i).dat_y(remap_idx(i).after,:);
    elseif config.normalization == 2
        fin_res(i).dat  = sep_results(i).dat(remap_idx(i).after,:);
    end
    fin_res(i).resp     = sep_results(i).resp(remap_idx(i).after);
end

% put both handpositions together in the last fin_res entry..
fin_res(i+1).dat    = [fin_res(1).dat; fin_res(2).dat];
fin_res(i+1).resp   = [fin_res(1).resp fin_res(2).resp];
