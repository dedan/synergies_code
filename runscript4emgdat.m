function runscript4emgdat


config.data_dir = '/Volumes/LAB/vega/data/';
config.ses_data = '/Volumes/LAB/vega/VegaSession';
config.monk     = 'v';
config.eddir    = '/Volumes/LAB/vega/EDfiles/';
config.e2add    = 'e';
config.outpath  = '/Volumes/LAB/vega/EMGdat/EMG';



vdir = dir([config.data_dir config.monk '*']);
                
vdir= sortdirs( vdir);

for i=1:length(vdir),
    
    curdir = char(vdir(i));
    if ~exist(['EMGdat' filesep 'EMG' curdir '.mat'], 'file'),

        disp(curdir);
        [em2take,e1,e2] = findEMGchannels( curdir, config);
        if ~isempty(em2take),
            [chdata,emgpsth] = muscleSyn( curdir, e1,e2, em2take, config); %#ok<NASGU>
            if ~isempty(chdata),
                save(['EMG' curdir], 'chdata', 'emgpsth');
                disp(['save EMG' curdir ' chdata emgpsth']);
            else
                disp('No available data');
            end
        end
    else
        disp(['EMGdat' filesep 'EMG' curdir '.mat skipped because already exists']);
    end
end


function odir = sortdirs( idir )

Ndates  = NaN(length(idir),2);
odir    = cell(length(idir));
for i=1:length(idir),
    curname = char(idir(i).name);
    Nday = str2double( curname(2:3));
    Nmonth = str2double(curname(4:5));
    Ndates(i,1) = Nmonth;
    Ndates(i,2) = Nday;
    odir(i) = {curname};
end
[~,ilist]=sortrows(Ndates);
odir = (odir(ilist));
