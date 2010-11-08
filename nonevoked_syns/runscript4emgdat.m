function runscript4emgdat


% config.data_dir = '/Volumes/LAB/vega/data/';
% config.ses_data = '/Volumes/LAB/vega/VegaSession';
% config.monk     = 'v';
% config.eddir    = '/Volumes/LAB/vega/EDfiles/';
% config.e2add    = 'e';
% config.outpath  = '/Volumes/LAB/vega/EMGdat/EMG';


config.data_dir = 'F:\vega\data\';
config.ses_data = 'F:\vega\VegaSession1';
config.monk     = 'v';
config.eddir    = 'F:\vega\EDfiles\';
config.e2add    = 'e';
config.outpath  = 'F:\vega\EMGdat\EMG';



vdir = dir([config.data_dir config.monk '*']);
                
        
vdir= sortdirs( vdir);

for i=1:length(vdir),
    
    if i < 19
        continue
    end
    
    curdir = char(vdir(i));
    if ~exist(['EMGdat' filesep 'EMG' curdir '.mat'], 'file'),

        disp(curdir);
        [em2take,e1,e2] = findEMGchannels( curdir, config);
        if ~isempty(em2take),
            [chdata,emgpsth] = muscleSyn( curdir, e1,e2, em2take, config);
            if ~isempty(chdata),
                cmnd  = ['save ' config.outpath curdir ' chdata emgpsth'];
                disp(cmnd);
                eval(cmnd);
            else
                disp('No available data');
            end
        end
    end
end


function odir = sortdirs( idir )

for i=1:length(idir),
    curname = char(idir(i).name);
    Nday = str2num( curname(2:3));
    Nmonth = str2num(curname(4:5));
    Ndates(i,1) = Nmonth;
    Ndates(i,2) = Nday;
    odir(i) = {curname};
end
[indx,ilist]=sortrows(Ndates);
odir = (odir(ilist));
