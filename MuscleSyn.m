function [chdata, emgpsth] = MuscleSyn( sessname, f1, f2, chnls, config)

chdata = [];
emgpsth = struct;

load([config.path config.monk 'session']);
emgdir = [config.path 'data' filesep sessname filesep 'mat' filesep ];
eddir  = [config.path 'EDfiles' filesep];


indx        = 1;
point2empty = [];
edfiles     = cell(0);

emgname = cell(1,length(f1:f2));
for i=f1:f2,
    edname = extract_edname( config.monk(1), subss, sessname,  i, config.e2add );
    if ~isempty(edname),
        if i< 10,
            str2add  = ['00' num2str(i)];
        elseif i< 100,
            str2add  = ['0' num2str(i)];
        end
        emgname(i-f1+1) = {[emgdir sessname str2add '_emg.mat']};
        if ~isempty(edname),
            edfiles(i-f1+1)     = {[eddir edname]};
        else
            point2empty(indx)   = i-f1+1; %#ok<AGROW>
            indx                = indx+1;
        end
    end
end

if isempty(edfiles),
    chdata  = [];
    emgpsth = [];
    return;
end

if ~isempty(point2empty),
    for j=1:length(point2empty),
        curpos = point2empty(j);
        disp(['For emg files: ' char(emgname(curpos)) '--> no edfiles was found']);
    end
    indx = ones(size(f1:f2));
    indx(point2empty) = 0;
    indx = find(indx);
    edfiles = edfiles(indx);
    emgname = emgname(indx);
end

Data = loadEMGdata( emgname, edfiles, chnls, config);
Nch  = length(Data.channel);
if Nch == 0,
    disp('No emg data was found');
    return;
end

trg = Data.channel(1).Target;               % target of trial
Ntr = unique(trg);
Ntr = Ntr(Ntr > 0);                         % list of targets
hnd = Data.channel(1).hand_position;        % hand position is the same for all channels
Nh  = unique(hnd(trg > 0));                 % list of hand positions

chdata = struct;
for k=1:length(Nh),
    
    % now computing the mean bckground level per channel
    for i=1:Nch,
        emgpsth(i).hand(k).bck_amp = mean(Data.channel(i).bck_amp);
    end
    
    for j=1:length(Ntr),
        jindx = find(trg == Ntr(j) & hnd == Nh(k));
        chdata(k).N(j) = length(jindx);
        for i=1:Nch,
            if length(jindx) >= 1,
                
                % this is the relevant data for further investigation
                chdata(k).amp{i,j}     = Data.channel(i).amp(jindx);
                chdata(k).bck_amp{i,j} = Data.channel(i).bck_amp(jindx);
                
                if length(jindx)> 1,
                    emgpsth(i).hand(k).target(j,:) = mean(Data.channel(i).signal(jindx,:));
                else
                    emgpsth(i).hand(k).target(j,:) = (Data.channel(i).signal(jindx,:));
                end
            else
                chdata(k).mat(i,j) = 0;
            end
        end;
    end
end



function   edname = extract_edname( monk,  subss, sessname, findx , e2add)

i = 1;
edname = '';
while i < length(subss),
    if strcmp( char(subss(i).Session), sessname ), % this is the right session,
        if ismember( findx, subss(i).Files(1):subss(i).Files(end)),
            ID = subss(i).ID;
            subID = subss(i).SubSess;
            ednum = find(  subss(i).Files(1):subss(i).Files(end) == findx);
            if ID < 10,
                str1 = ['0' num2str(ID)];
            else
                str1 = num2str(ID);
            end
            if subID < 10,
                str2 = ['0' num2str(subID)];
            else
                str2 = num2str(subID);
            end
            edname = [monk str1 str2 'e' e2add '.' num2str(ednum) '.mat'];
            return;
        end
    end
    i = i+1;
end
disp( 'EDname not found' );



