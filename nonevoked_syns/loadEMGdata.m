function tuning = loadEMGdata( emgfiles, edfiles, channels )
tuning = [];
bhv = 'TO';
bck_bhv = 'SC';
pre = -500;
post = 1000;
win = [0 500];% time relative to torque onset TO
total_trials = 0;
bck_win = [-500 0]; % time relative to spatial cue (SC)

for i=1:length(emgfiles),
    if isempty(emgfiles{i}),
        curfile = '';
        curbhv = '';
    else
        curfile = char(emgfiles(i));
        curbhv = char(edfiles(i));
    end
    try,
        eval(['tmp = load(curfile,''' 'EMG' num2str(channels(1)) '_KHz' '''' ');']);
    catch,
        if exist(curfile,'file'),
            disp([curfile ' Corrupted!!!']);
            curfile = '';
        end
    end
    if exist( curfile,'file') & exist(curbhv, 'file')
        bhvdata  = load(curbhv);
        disp(curbhv);
        try
            hand_position=bhvdata.hand_position;
            handFlag=1;
        catch
            warning('No hand position variable found')
            handFlag=0;
        end
        if ~isfield( bhvdata, 'bhvStat'),
            trials=locate_trials(bhvdata,'all');
            for tri=1:size(trials,1),
                [direction, TargetNum]=getDir(bhvdata,trials(tri,:));
                targets(tri) = TargetNum;
            end
        else,
            targets = [];
            bhvStat = bhvdata.bhvStat;
            if ~isempty(bhvStat),
                % STEPHAN: here we select trials based on some behavioral
                % criteria of reaction time, {-200-500) Movement time
                % (500-1500) and some angular deviation (<35). I think we
                % could relax some of these criteria) especially the last
                % one
                itrials=(find(bhvStat(:,1)>=-200 & bhvStat(:,1)<=500 & bhvStat(:,2)<=1500  & bhvStat(:,2)>=500 & abs(bhvStat(:,5))<=35));
                trials = bhvdata.trials(itrials,1:2);
                if isfield(bhvdata, 'targets'),
                    targets = bhvdata.targets(itrials);
                    if ~isempty(find(targets==0)),
                        disp('here');
                    end
                else
                    disp('creating target data');
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
        end
        trialCounter = 0;
        if ~isempty(trials),
            for m=1:length(channels)
                chstr=num2str(channels(m));
                eval(['load(curfile,''' 'EMG' chstr '''' ')']);
                eval(['load(curfile,''' 'EMG' chstr '_KHz' '''' ')']);
                eval(['EMG=EMG' chstr ';']);
                eval(['Fs = EMG' chstr '_KHz;']);
                EMG = abs(EMG-mean(EMG));
                df = 100;
                EMG = decimate( EMG, df);
                Fs = Fs/df;
                Lemg = length(EMG);
                time=[1:Lemg]/Fs; %%%change this 5!!!!!!!!!!!!!!!!\
                for k=1:size(trials,1)
                    trialCounter=trialCounter+1;
                    TargetNum = targets(k);
                    direction = TargetNum;
                    %                     [direction, TargetNum]=getDir(bhvdata,trials(k,:));
                    refTime=getRefTime(bhvdata,trials(k,:),bhv, TargetNum);
                    %neuronOffset=neurons(1)-1;
                    index=find(time>=refTime+pre & time<=refTime+post);
                    index2=find(time>=refTime+win(1) & time<=refTime+win(2));
                    if mod(length(index),2)~=0
                        index=index(1:end-1);
                    end
                    signal = EMG(index);
                    signal2 = EMG(index2);
                    amp=mean(abs(signal2));
                    tuning.channel(m).dir(trialCounter+total_trials)=direction;
                    tuning.channel(m).Target(trialCounter+total_trials)=TargetNum;
                    tuning.channel(m).amp(trialCounter+total_trials)=amp;
                    tuning.channel(m).signal(trialCounter+total_trials,1:length(signal))=abs(signal);
                    
                    % we are now computing the backcground emg level - 500
                    % ms before cue onset
                    
                    refTime=getRefTime(bhvdata,trials(k,:),bck_bhv, TargetNum);
                    index=find(time>=refTime+bck_win(1) & time<=refTime+bck_win(2));
                    bck_signal = EMG(index);
                    bck_amp=mean(abs(signal));
                    tuning.channel(m).bck_amp(trialCounter+total_trials) = bck_amp;
            
                    % %                     if length(signal) < size(tuning.channel(m).signal,2),
                    % %                         disp('here');
                    % %                     end
                    if handFlag
                        tuning.channel(m).hand_position(trialCounter+total_trials)=hand_position;
                    end
                end
            end
            eval(['clear EMG' chstr ' EMG'])
        end
    else
        trialCounter = 0;
    end
	total_trials = total_trials + trialCounter;
end

        
            
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trials=locate_trials(fileDat,mod)

switch mod
    case 'all'
        %strt_index=find(ismember(fileDat.events_code,hex2dec(['100';'200';'201';'206';'207';'209'])));        
        strt_index=find(ismember(fileDat.events_code,hex2dec(['100';'209']))); 
        for i=1:length(strt_index)
            end_index(i)=min(find(ismember(fileDat.events_code(strt_index(i):end),hex2dec('101'))))+strt_index(i)-1;
            if isempty(end_index)
                strt_index=strt_index(1:end-1);
            end
        end
    case '100'
        strt_index=find(ismember(fileDat.events_code,hex2dec('100')));
        for i=1:length(strt_index)
            end_index(i)=min(find(ismember(fileDat.events_code(strt_index(i):end),hex2dec('101'))))+strt_index(i)-1;
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
targetCode=trialEvents(find(ismember(trialEvents,hex2dec(['20';'21';'22';'23';'24';'25';'26';'27']))));
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

i=find(ismember(events,hex2dec(['20';'21';'22';'23';'24';'25';'26';'27'])));
targetNum=events(i)-31;