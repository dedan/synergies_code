function [chdata, emgpsth] = MuscleSyn( sessname, config)

emgpsth = struct;
chdata  = [];

emgdir = [config.path 'data' filesep sessname filesep 'mat' filesep ];
eddir  = [config.path 'EDfiles' filesep];


% load the info file
load([config.path 'data' filesep sessname filesep 'info' filesep sessname '_param.mat']);

n_emg_files = SESSparam.SubSess(end).Files(2);
emgname     = cell(1,n_emg_files);
edfiles     = cell(1,n_emg_files);


for i = 1:n_emg_files
    
    % construct the name of the EDfile
    str1  = sprintf('%02d', DDFparam.ID);
    str2  = '';
    found = false;
    j     = 1;
    
    while ~found && j <= length(SESSparam.SubSess)
        ednum = find(SESSparam.SubSess(j).Files(1):SESSparam.SubSess(j).Files(2) == i);
        if ~isempty(ednum)
            str2  = sprintf('%02d', j);
            found = true;
        end
        j = j + 1;
    end
    
    % NOTE this can be removed if I never find it in the log
    if isempty(str2)
        disp('EDname not found');
    end
    
    edname = [sessname(1) str1 str2 'e' config.e2add '.' num2str(ednum) '.mat'];
    
    str2add = sprintf('%03d', i);
    emgname(i) = {[emgdir sessname str2add '_emg.mat']};
    edfiles(i) = {[eddir edname]};
end

% NOTE this can be removed if I never find it in the log
if isempty(edfiles)
    disp('edfiles is empty');
end


Data = loadEMGdata( emgname, edfiles, config);
if isempty(Data)
    disp('problem in loadEMG');
    return;
end

Nch  = length(Data.channel);
if Nch == 0,
    disp('No emg data was found');
    return;
end

trg       = Data.channel(1).Target;               % target of trial
Ntr       = unique(trg);
Ntr       = Ntr(Ntr > 0);                         % list of targets
hnd       = Data.channel(1).hand_position;        % hand position is the same for all channels
Nh        = unique(hnd(trg > 0));                 % list of hand positions

for k=1:length(Nh),
    
    chdata(k).channels = Data.channels; %#ok<AGROW>
    chdata(k).id = DDFparam.ID; %#ok<AGROW>
    
    % now computing the mean bckground level per channel
    for i=1:Nch,
        emgpsth(i).hand(k).bck_amp = mean(Data.channel(i).bck_amp);
    end
    
    for j=1:length(Ntr),
        jindx = find(trg == Ntr(j) & hnd == Nh(k));
        
        tmp_amp = NaN(Nch, length(jindx));
        tmp_bck = NaN(Nch, length(jindx));
        for i=1:Nch,
            if length(jindx) >= 1,
                
                % this is the relevant data for further investigation
                tmp_amp(i,:) = Data.channel(i).amp(jindx);
                tmp_bck(i,:) = Data.channel(i).bck_amp(jindx);
                
                if length(jindx)> 1,
                    emgpsth(i).hand(k).target(j,:) = mean(Data.channel(i).signal(jindx,:));
                else
                    emgpsth(i).hand(k).target(j,:) = (Data.channel(i).signal(jindx,:));
                end
            end
        end;
        chdata(k).amp{j}     = tmp_amp; %#ok<AGROW>
        chdata(k).bck_amp{j} = tmp_bck; %#ok<AGROW>
    end
end
