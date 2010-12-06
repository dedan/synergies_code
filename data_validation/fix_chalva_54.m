
% chalva sessions 54 the connector was switched. 
% channels 1-8 were switched with 9-16
function fix_chalva_54(path)

path = [path 'chalva' filesep 'EMGdat' filesep 'EMGc040810.mat'];
load(path);

idx  = [9:16 1:8];

emgpsth   = emgpsth(idx); %#ok<NODEF,NASGU>
chdata.pd = chdata.pd(idx); %#ok<NODEF>
chdata.p1 = chdata.p1(idx);
chdata.p2 = chdata.p2(idx);

for i = 1:length(chdata.amp)
    chdata.amp{i} = chdata.amp{i}(idx,:);
end

save(path, 'chdata', 'emgpsth');