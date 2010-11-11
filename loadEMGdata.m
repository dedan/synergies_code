function tuning = loadEMGdata( emgfiles, edfiles, channels, config )

df      = 100;
tuning  = [];
bhv     = 'TO';
bck_bhv = 'SC';
pre     = -500;
post    = 1000;
win     = [0 500];  % time relative to torque onset TO
bck_win = [-500 0]; % time relative to spatial cue (SC)
total_trials = 0;


for i=1:length(emgfiles)

    curfile = '';
    curbhv  = '';
    if ~isempty(emgfiles{i})
        curfile = char(emgfiles(i));
        curbhv  = char(edfiles(i));
    end
    
    % check whether file contains channel data
    try
        load(curfile,['EMG' num2str(channels(1)) '_KHz']);
    catch err
        if config.verbose
            disp(err);
        end
        if exist(curfile,'file'),
            disp([curfile ' Corrupted!!!']);
            curfile = '';
        end
    end
    
    if exist(curfile,'file') && exist(curbhv, 'file')
        bhvdata = load(curbhv);
        if config.verbose
            disp(curbhv);
        end
        try
            hand_position   = bhvdata.hand_position;
            handFlag        = 1;
        catch %#ok<CTCH>
            warning('muscle_syn:no_hand', 'No hand position variable found')
            handFlag = 0;
        end
        
        % compute targets if bhvStat does not exist
        if ~isfield( bhvdata, 'bhvStat')
            if 1    % config.verbose
                disp 'computing targets' %just want to check whether this condition ever occurs
            end
            trials  = locate_trials(bhvdata,'all');
            targets = NaN(1,size(trials,1));
            for tri=1:size(trials,1),
                [~, TargetNum] = getDir(bhvdata,trials(tri,:));
                targets(tri)   = TargetNum;
            end
        elseif ~isempty(bhvdata.bhvStat)
            bhvStat = bhvdata.bhvStat;
            % STEPHAN: here we select trials based on some behavioral
            % criteria of reaction time, {-200-500) Movement time
            % (500-1500) and some angular deviation (<35). I think we
            % could relax some of these criteria) especially the last
            % one
            itrials =(find(bhvStat(:,1)>=-200 & bhvStat(:,1)<=500 & bhvStat(:,2)<=1500  & bhvStat(:,2)>=500 & abs(bhvStat(:,5))<=35));
            trials  = bhvdata.trials(itrials,1:2);
            if isfield(bhvdata, 'targets'),
                if config.verbose
                    disp('found target data');
                end
                targets = bhvdata.targets(itrials);
            else
                if config.verbose
                    disp('creating target data');
                end
                targets = NaN(1,size(trials,1));
                for ie = 1:size(trials,1),
                    tr1 = trials(ie,1);
                    tr2 = trials(ie,2);
                    events=bhvdata.events_code(tr1:tr2);
                    targets(ie)=returnTarget(events);
                end
            end
        else
            trials = [];
        end
        
        if ~isempty(trials),
            for m=1:length(channels)
                
                % load the data
                chstr   = num2str(channels(m));
                EMG     = load(curfile, ['EMG' chstr]);
                EMG     = EMG.(['EMG' chstr]);
                Fs      = load(curfile, ['EMG' chstr '_KHz']);
                Fs      = Fs.(['EMG' chstr '_KHz']);
                
                % preprocessing
                EMG     = abs(EMG-mean(EMG));
                EMG     = decimate( EMG, df);
                Fs      = Fs/df;
                Lemg    = length(EMG);
                time    = (1:Lemg)/Fs; 
                
                for k=1:size(trials,1)
                    direction       = targets(k);
                    refTime         = getRefTime(bhvdata,trials(k,:),bhv, direction);
                    
                    % select signals in the windows
                    index           = find(time>=refTime+pre & time<=refTime+post);
                    index2          = time>=refTime+win(1) & time<=refTime+win(2);
                    if mod(length(index),2)~=0
                        index=index(1:end-1);
                    end
                    signal  = EMG(index);
                    signal2 = EMG(index2);
                    tuning.channel(m).dir(k+total_trials)    = direction;
                    tuning.channel(m).Target(k+total_trials) = direction;
                    tuning.channel(m).amp(k+total_trials)    = mean(abs(signal2));
                    tuning.channel(m).signal(k+total_trials,1:length(signal)) = abs(signal);
                    
                    % we are now computing the backcground emg level - 500 ms before cue onset
                    refTime     = getRefTime(bhvdata,trials(k,:),bck_bhv, direction);
                    index       = time>=refTime+bck_win(1) & time<=refTime+bck_win(2);
                    bck_signal  = EMG(index);
                    bck_amp     = mean(abs(bck_signal));
                    tuning.channel(m).bck_amp(k+total_trials) = bck_amp;
            
                    if handFlag
                        tuning.channel(m).hand_position(k+total_trials) = hand_position;
                    end
                end
            end
            eval(['clear EMG' chstr ' EMG'])
        end
    end
	total_trials = total_trials + size(trials,1);
end

        
            
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trials=locate_trials(fileDat,mod)


switch mod
    case 'all'      
        strt_index=find(ismember(fileDat.events_code,hex2dec(['100';'209']))); 
        end_index = NaN(1,length(strt_index));
        for i=1:length(strt_index)
            end_index(i)=find(ismember(fileDat.events_code(strt_index(i):end),hex2dec('101')), 1 )+strt_index(i)-1;
            if isempty(end_index)
                strt_index=strt_index(1:end-1);
            end
        end
    case '100'
        strt_index=find(ismember(fileDat.events_code,hex2dec('100')));
        end_index = NaN(1,length(strt_index));
        for i=1:length(strt_index)
            end_index(i)=find(ismember(fileDat.events_code(strt_index(i):end),hex2dec('101')), 1 )+strt_index(i)-1;
            if isempty(end_index)
                strt_index=strt_index(1:end-1);
            end
        end
end
if isempty(strt_index),
    disp('no trials were found');
    trials = [];
    return;
end
trials=[strt_index end_index'];


function [rad,TargetNum]=getDir(fileDat,trial)

trialEvents=fileDat.events_code(trial(1):trial(2));
targetCode=trialEvents(ismember(trialEvents,hex2dec(['20';'21';'22';'23';'24';'25';'26';'27'])));
switch targetCode
    case 32
        rad=pi/2;
        TargetNum=1;
    case 33
        rad=pi/4;
        TargetNum=2;
    case 34
        rad=0;
        TargetNum=3;
    case 35
        rad=-pi/4;
        TargetNum=4;
    case 36
        rad=-pi/2;
        TargetNum=5;
    case 37
        rad=-3*pi/4;
        TargetNum=6;
    case 38
        rad=pi;
        TargetNum=7;
    case 39
        rad=3*pi/4;
        TargetNum=8;
end

function refTime=getRefTime(fileDat,trial,bhv, TargetNum)

trialEvents=fileDat.events_code(trial(1):trial(2));
switch bhv
    case 'ST'
        indx=trial(1);
    case 'SC'
        code=hex2dec(['2' num2str(TargetNum-1)]);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'GO'
        code=hex2dec(['3' num2str(TargetNum-1)]);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'TO'
        code=hex2dec([num2str(TargetNum) 'A']);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'HO'
        code=hex2dec([num2str(TargetNum) 'B']);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'B2C'
        code=64;
        indx=find(trialEvents==code)+trial(1)-1;
    case 'HOF'
        code=hex2dec([num2str(TargetNum) 'C']);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'TOF'
        code=hex2dec([num2str(TargetNum) 'D']);
        indx=find(trialEvents==code)+trial(1)-1;
    case 'ET'
        indx=trial(2);
end
refTime=fileDat.events_time(indx);
 

function targetNum=returnTarget(events)

i= ismember(events,hex2dec(['20';'21';'22';'23';'24';'25';'26';'27']));
targetNum=events(i)-31;