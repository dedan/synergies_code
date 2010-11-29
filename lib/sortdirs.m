function odir = sortdirs( idir )

Ndates  = NaN(length(idir),2);
odir    = cell(length(idir));
for i=1:length(idir),
    curname     = char(idir(i).name);
    Nday        = str2double( curname(2:3));
    Nmonth      = str2double(curname(4:5));
    Ndates(i,1) = Nmonth;
    Ndates(i,2) = Nday;
    odir(i)     = {curname};
end
[~,ilist]   =sortrows(Ndates);
odir        = (odir(ilist));