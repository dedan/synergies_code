function data = load_emg(edfiles, emgfiles, config )

total_trials = 0;
data         = [];
all          = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];


c2take   = get_channels2take(emgfiles{1});
channels = zeros(1,length(all));
for j = 1:length(c2take)
    channels(all == c2take(j)) = 1;
end
data.channels = channels;


for i = 1:length(emgfiles)
    bhvdata = load(edfiles{i});

    % select trials with certain properties
    if isempty(bhvdata.bhvStat)
        trials = [];
    else
        itrials = bhvdata.bhvStat(:,1) >= config.t_react(1) & ...
            bhvdata.bhvStat(:,1) <= config.t_react(2) & ...
            bhvdata.bhvStat(:,2) >= config.t_move(1) &  ...
            bhvdata.bhvStat(:,2) <= config.t_move(2) &  ...
            abs(bhvdata.bhvStat(:,5)) <= config.ang_div;
        trials  = bhvdata.trials(itrials,1:2);
    end
    
    if isempty(trials)
        continue;
    end        
    
    for j = 1:length(c2take)
        
        % load the data
        chstr   = num2str(c2take(j));
        EMG     = load(emgfiles{i}, ['EMG' chstr]);
        EMG     = EMG.(['EMG' chstr]);
        Fs      = load(emgfiles{i}, ['EMG' chstr '_KHz']);
        Fs      = Fs.(['EMG' chstr '_KHz']);
        
        % preprocessing
        EMG     = abs(EMG-mean(EMG));
        EMG     = decimate( EMG, config.df);
        Fs      = Fs/config.df;
        Lemg    = length(EMG);
        time    = (1:Lemg)/Fs;
    
        for k = 1:size(trials,1)
            [dir, target] = getDir(bhvdata, trials(k,:));
            refTime       = getRefTime(bhvdata,trials(k,:),config.bhv, target);
            
            if refTime > max(time)
                disp('refTime Problem');
                data = [];
                return;
            end
            
            % select signals in the windows
            index           = find(time>=refTime+config.pre & time<=refTime+config.post);
            index2          = time>=refTime+config.win(1) & time<=refTime+config.win(2);
            if mod(length(index),2)~=0
                index=index(1:end-1);
            end
            signal  = EMG(index);
            signal2 = EMG(index2);
            data.channel(j).dir(k+total_trials)    = dir;
            data.channel(j).Target(k+total_trials) = target;
            data.channel(j).amp(k+total_trials)    = mean(abs(signal2));
            data.channel(j).signal(k+total_trials,1:length(signal)) = abs(signal);
            
            % we are now computing the backcground emg level - 500 ms before cue onset
            refTime     = getRefTime(bhvdata,trials(k,:),config.bck_bhv, target);
            index       = time>=refTime+config.bck_win(1) & time<=refTime+config.bck_win(2);
            bck_signal  = EMG(index);
            bck_amp     = mean(abs(bck_signal));
            data.channel(j).bck_amp(k+total_trials) = bck_amp;
            
            if isfield(bhvdata, 'hand_position')
                data.channel(j).hand_position(k+total_trials) = bhvdata.hand_position;
            else
                warning('muscle_syn:no_hand', 'No hand position variable found')
            end
        end
    end
    total_trials = total_trials + size(trials,1);
end


if total_trials == 0
    disp('all trials empty');
    data = [];
    return
end

u = unique(data.channel(i).hand_position);
u = u(u ~= 0);
for j = 1:length(u)
    data.pd = zeros(length(u), length(all));
    data.p1 = zeros(length(u), length(all));
    data.p2 = zeros(length(u), length(all));
    
    for i = 1:length(c2take)
        tmp                 = data;
        tmp.channel         = data.channel(i);
        indx                = find(data.channel(i).hand_position == u(j));
        tmp.channel.dir     = data.channel(i).dir(indx);
        tmp.channel.amp     = data.channel(i).amp(indx);
        tmp.channel.signal  = data.channel(i).signal(indx,:);
        tmp.channel.hand_position = data.channel(i).hand_position(indx);
        
        pd = get_pd( tmp);
        p1 = sig_dir_emg(tmp.channel, config);
        p2 = anova1(tmp.channel.amp, tmp.channel.dir,'off');
        
        data.pd(j,i)     = pd;
        data.p1(j,i)     = p1;
        data.p2(j,i)     = p2;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c2take = get_channels2take(emgfile)

vrs     = who('-file', emgfile);
indx    = 1;
for j=1:length(vrs),
    curvar = char(vrs(j));
    if ~isempty(findstr(curvar,'EMG')) && isempty(findstr(curvar,'KHz')),
        c2take(indx)  = sscanf(curvar,'EMG%d'); %#ok<AGROW>
        indx          = indx+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pd = get_pd(data )

u = unique(data.channel.dir);
x = NaN(1,length(u));
y = NaN(1,length(u));

for i=1:length(u);
    tmp  = data.channel.dir == u(i);
    x(i) = u(i);
    y(i) = mean(data.channel.amp(tmp));
end

cx  = cos(x);
sx  = sin(x);
mcx = sum(y.*cx)/sum(y);
msx = sum(y.*sx)/sum(y);
pd  = atan(msx/mcx)+(mcx<0)*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = sig_dir_emg(data, config)

teta=[90 45 0 -45 -90 -135 180 135]*pi/180;

if sum(data.amp) ==0,
    p=1;
    return;
end;

u           = unique(data.dir);
teta2take   = teta(ismember(teta, u));
randP       = round(rand(length(data.amp)*config.n_btstrp,1)*length(data.amp)+0.5); 
indx        = reshape(randP,[],config.n_btstrp);
bootRates   = data.amp(indx);

count = 1;
means = zeros(1,length(u));
 for i = 1:length(u)
     
     tmpRates   = data.amp(data.dir == u(i));
     means(i)   = mean(tmpRates);
     ntmp       = length(find(data.dir == u(i)));
     vec2take   = bootRates(count:count+ntmp-1,:);
     if size(vec2take,1) == 1,
         bM(i,:)=vec2take; %#ok<AGROW>
     else
         bM(i,:)=mean(vec2take); %#ok<AGROW>
     end;
     count=count+ntmp;
end

Rx      = cos(teta2take)*means';
Ry      = sin(teta2take)*means';
R       = sqrt(Rx.^2+Ry.^2);
Rx      = cos(teta2take)*bM;
Ry      = sin(teta2take)*bM;
Rboot   = sqrt(Rx.^2+Ry.^2);
p       = length(find(Rboot>R))/config.n_btstrp;
