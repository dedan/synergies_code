function [chdata, emgpsth] = get_emg4sess(sessname, config)

% initialize return values
emgpsth = struct;
chdata  = [];

% set the data folders and load info file
emgdir = [config.path 'data' filesep sessname filesep 'mat' filesep ];
eddir  = [config.path 'EDfiles' filesep];
%info = load([config.path 'data' filesep sessname filesep 'info' filesep sessname '_param.mat']);
info = load([config.path 'info_files' filesep sessname '_param.mat']);


[edfiles, emgfiles] = get_files2run(info, eddir, emgdir, sessname, config);

if length(edfiles) ~= length(emgfiles)
    disp('some files must be missing, stop working on this session');
    return;
end

if isempty(emgfiles)
    disp('no valid subsessions found, stop working this session');
    return
end


data = load_emg(edfiles, emgfiles, config);

if isempty(data)
    disp('data empty, stop working this session');
    return
end

take = find(data.channels);
trg  = data.channel(take(1)).Target;         % target of trial
Ntr  = unique(trg);
Ntr  = Ntr(Ntr > 0);                         % list of targets
hnd  = data.channel(take(1)).hand_position;  % hand position is the same for all channels
Nh   = unique(hnd(trg > 0));                 % list of hand positions
Nh   = Nh(Nh ~= 0);

for k=1:length(Nh),
    
    chdata(k).channels  = logical(data.channels);  %#ok<*AGROW>
    chdata(k).id        = info.DDFparam.ID; 
    chdata(k).pd        = data.pd(k,:);
    chdata(k).p1        = data.p1(k,:);
    chdata(k).p2        = data.p2(k,:);
    chdata(k).trials    = data.trials;

    
    % now computing the mean bckground level per channel
    for i=find(data.channels)
        emgpsth(i).hand(k).bck_amp = mean(data.channel(i).bck_amp);
    end
    
    for j=1:length(Ntr),
        jindx = find(trg == Ntr(j) & hnd == Nh(k));
        
        tmp_amp = NaN(length(find(data.channels)), length(jindx));
        tmp_bck = NaN(length(find(data.channels)), length(jindx));
        
        c_inds = find(data.channels); 
        for i = 1:length(c_inds)
            if length(jindx) >= 1,
                
                % this is the relevant data for further investigation
                tmp_amp(i,:) = data.channel(c_inds(i)).amp(jindx);
                tmp_bck(i,:) = data.channel(c_inds(i)).bck_amp(jindx);
                
                if length(jindx)> 1,
                    emgpsth(c_inds(i)).hand(k).target(j,:) = mean(data.channel(c_inds(i)).signal(jindx,:));
                else
                    emgpsth(c_inds(i)).hand(k).target(j,:) = (data.channel(c_inds(i)).signal(jindx,:));
                end
            end
        end;
        
        % give all data the same size (as if recorded from all channels)
        chdata(k).amp{j}     = zeros(length(chdata(k).channels), size(tmp_amp,2));
        chdata(k).bck_amp{j} = zeros(length(chdata(k).channels), size(tmp_bck,2));
        chdata(k).amp{j}(chdata(k).channels,:)      = tmp_amp;
        chdata(k).bck_amp{j}(chdata(k).channels,:)  = tmp_bck;
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [edfiles, emgfiles] = get_files2run(info, eddir, emgdir, curname, config)

edfiles  = cell(0);
emgfiles = cell(0);
for i=1:length(info.SESSparam.SubSess),
    subs_files = info.SESSparam.SubSess(i).Files(1):info.SESSparam.SubSess(i).Files(end);
    stimflag   = test4stim(info.SESSparam.SubSess(i));
    configflag = test4other(info.SESSparam.fileConfig, subs_files);
    if config.also_with_stim || (~stimflag && ~configflag)
        new      = get_ed_files(eddir, i, info.DDFparam.ID, subs_files, curname, config);
        edfiles  = horzcat(edfiles, new); %#ok<AGROW>
        new      = get_emg_files(emgdir, curname, subs_files );
        emgfiles = horzcat(emgfiles, new); %#ok<AGROW>
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimflag = test4stim(subs)

stimflag = 0;
if isfield(subs, 'Electrode')
    for i=1:length(subs.Electrode)
        if isfield(subs.Electrode(i).Stim, 'Flag') && subs.Electrode(i).Stim.Flag
            stimflag = 1;
            return;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function configflag = test4other(fileConfig, files)

configflag = 0;

for i=files,
    if (isfield(fileConfig, 'SCPStim')  && fileConfig(i).SCPStim) || ...
       (isfield(fileConfig, 'PRB')      && fileConfig(i).PRB) || ...
       (isfield(fileConfig, 'DBS')      && fileConfig(i).DBS)
        configflag = 1;
        return
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  edfiles = get_ed_files(eddir, subindx, sessID, files, curname, config )

sess_str = sprintf('%02d', sessID);
sub_str  = sprintf('%02d', subindx);

edfiles = cell(1,length(files));
for i=1:length(files),
    name2add = [eddir curname(1) sess_str sub_str config.e2add '.' num2str(i) '.mat'];
    if exist(name2add,'file'),
        edfiles{i} = name2add;
    else
        disp(['missing edfile: ' name2add ]);
    end
end
edfiles = edfiles(~cellfun('isempty',edfiles));
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emgfiles = get_emg_files(matdir, curname, files )

emgfiles = cell(1,length(files));
for i=1:length(files)
    ex = sprintf('%03d', files(i)); 
    name2add = [matdir curname ex '_emg.mat'];
    if exist(name2add,'file'),
        emgfiles{i} = name2add;
    else
        disp(['missing emgfile: ' name2add ]);
    end
end
emgfiles = emgfiles(~cellfun('isempty',emgfiles));

