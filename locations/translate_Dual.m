function [dx,dy] = translate_Dual( iel, curdir, x,y);
% some translation with respect to date of recording session!??


DD = str2num(curdir(2:3));
MM = str2num(curdir(4:5));

if (DD >= 15 & MM ==2) | (MM > 2),
    if iel == 1, % left electrode
        dx = x + 6.5;
        dy = y +1;
    elseif iel == 2,
        dx =x-4;
        dy = y-2;
    end
elseif DD==14 & MM==2,
    if iel == 1, % left electrode
        dx = x + 7.5;
        dy = y +.5;
    elseif iel == 2,
        dx =x-5;
        dy = y-1;
    end
elseif DD>= 7 & MM== 2,
        if iel == 1, % left electrode
        dx = x + 10.8;
        dy = y -2;
    elseif iel == 2,
        dx =x-11.5;
        dy = y+0.5;
    end
else
    disp('Stop!!');
end
