% here we add the mapping info to the response structure. The mapping
% include a 1x2 cell array. The components are:
% 1. the site according to the stimulation threshold (<= 10 m1 >= 15 pm and
% between ambiguous).
% 2. the site according to the anatomical location (relative to sulcus as
% anatomical landmarks

function add_mapinfo(indir, monks)

    load(['data' filesep 'mismatches.mat']);
    dat = r(2:end,:);
    for i=1:size(dat,1),
        fnm = lower(dat{i,1});
        monkcode(i) = double(fnm(1));
    end

    for i=1:length(monks),
        curmonk = monks{i};
        load([indir 'evoked_data_' curmonk]);
        indx2take = find(monkcode == double(curmonk(1)));
        resps = addmap2resps(resps, dat(indx2take,:));
        save( [indir 'evoked_data_map_' curmonk], 'resps');
    end
    disp('m1 pm information added');
end



function resps = addmap2resps(resps, msmtchs)

    for i=1:length(resps),
        if isempty( resps(i).location.thr ),
            resps(i).location.thr = nan;
        end
        if ischar(resps(i).location.thr),
            resps(i).location.thr = str2double( resps(i).location.thr);
        end
        mindx =  find_incell( msmtchs,2, resps(i).id, 3, resps(i).electrode  );
        if mindx,
            resps(i).mapsite = get_mismatch_data( msmtchs(mindx,:));
        else
            if resps(i).location.thr <= 10,
                resps(i).mapsite = {'m1' 'm1'};
            elseif resps(i).location.thr >= 15,
                resps(i).mapsite = {'pm' 'pm'};
            else
                disp(['MBG>>' char(resps(i).session) ',id=' num2str( resps(i).id) ',el=' num2str(resps(i).electrode) ', thr=' num2str(resps(i).location.thr)]);
            end
        end

    end
end


function indx =  find_incell( celldat,col1, n1, col2, n2 )

    num =  cell2mat(celldat(:,[col1 col2]));
    indx = find(num(:,1) == n1 & num(:,2) == n2);
end

function mapdat = get_mismatch_data( msdat )

    mapdat = msdat(end-1:end);
end



