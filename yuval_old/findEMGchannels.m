function [chns, e1, e2] = findEMGchannels( sessname, config )

load(config.ses_data);
emgdir = [config.data_dir sessname filesep 'mat' filesep ];


femgs = dir([emgdir '*_emg.mat']);
if isempty(femgs),
    disp('No EMG mat files');
    chns = [];
    e1 = 0;
    e2=0;
    return;
end

e1= 1;
e2 = length(femgs);

curname = [emgdir char(femgs(1).name)];
vrs = who('-file', curname);
indx =1;
for i=1:length(vrs),
    curvar = char(vrs(i));
    if ~isempty(findstr(curvar,'EMG')) & isempty(findstr(curvar,'KHz')),
        chns(indx) = sscanf(curvar,'EMG%d');
        indx = indx+1;
    end
end
