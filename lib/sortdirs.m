function dlist = sortdirsyp( alldirs )

indx = 1;
for i=1:length( alldirs  )
    if (~alldirs(i).isdir || ( strcmp(alldirs(i).name(1), '.'))), continue; end; % taking only directories
    curdir = char(alldirs(i).name);
    dlist(indx) = {curdir};
    tmp(indx,1) = str2num(curdir(6:7));
    tmp(indx,2) = str2num(curdir(4:5));
    tmp(indx,3) = str2num(curdir(2:3));
    indx = indx+1;

end
[tmp,i]=sortrows(tmp);
dlist = dlist(i);

