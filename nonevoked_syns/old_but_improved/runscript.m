function runscript


vdir = dir('/Volumes/LAB/vega/data/v*');
vdir= sortdirs( vdir);

for i=1:length(vdir),
    curdir = char(vdir(i));
    disp(curdir);
    [em2take,e1,e2] = findEMGchannels( curdir);
    if ~isempty(em2take),
        [chdata,emgpsth] = muscleSyn( curdir, e1,e2, em2take);
        if ~isempty(chdata),
            cmnd  = ['save EMG' curdir ' chdata emgpsth'];
            disp(cmnd);
            eval(cmnd);
        else
            disp('No available data');
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
