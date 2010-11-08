function [chdata, emgpsth] = MuscleSyn( sessname, f1, f2, chnls)

chdata = [];
PLOTFLAG = 1;

%NOTE pathes
load /Volumes/LAB/vega/VegaSession

emgdir = ['/Volumes/LAB/vega/data/' sessname '/mat/' ];
eddir = '/Volumes/LAB/vega/edfiles/';

indx = 1;
point2empty = [];
edfiles = cell(0);
emgname = cell(0);

for file_number = f1:f2
    edname = extract_edname(subss, sessname, file_number);
    if ~isempty(edname)
        str2add = [lead_zeros(file_number,3) num2str(file_number)];
        emgname(indx) = {[emgdir sessname str2add '_emg.mat']};
        edfiles(indx) = {[eddir edname]};
        indx = indx+1;
    else
        disp(['For emg files: ' emgdir sessname str2add '_emg.mat --> no edfiles was found']);
    end
end

if isempty(edfiles),
    chdata = [];
    emgpsth = [];
    return;
end


Data = loadEMGdata( emgname, edfiles, chnls);
Nch = length(Data.channel);
if Nch == 0,
    disp('No emg data was found');
    return;
end
trg = Data.channel(1).Target;
Ntr = unique(trg);
Ntr  = Ntr(find(Ntr > 0));
hnd = Data.channel(1).hand_position;
Nh = unique(hnd(find(trg > 0)));
for k=1:length(Nh),
    for j=1:length(Ntr),
        jindx = find(trg == Ntr(j) & hnd == Nh(k));
        chdata(k).N(j) = length(jindx);
        for i=1:Nch,
            if length(jindx) >= 1,
                chdata(k).mat(i,j) = mean(Data.channel(i).amp(jindx));
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
    %     subplot(4,4,i);
    %     bar(chdata(i,:));
end



function edname = extract_edname(subss, sessname, file_number)

i = 1;
edname = '';
while i < length(subss)
    if strcmp( char(subss(i).Session), sessname ) % this is the right session,
        if ismember( file_number, subss(i).Files(1):subss(i).Files(end))
            ID = subss(i).ID;
            subID = subss(i).SubSess;
            ednum = find(  subss(i).Files(1):subss(i).Files(end) == file_number);
            str1 = [lead_zeros(ID,2) num2str(ID)];
            str2 = [lead_zeros(subID,2) num2str(subID)];
            edname = ['v' str1 str2 'ee.' num2str(ednum) '.mat'];
            return;
        end
    end
    i = i+1;
end
disp( 'EDname not found' );



