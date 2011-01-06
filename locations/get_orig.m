% get the coordinates of an origin for a quadrant
function [x, y] = get_orig(qd, PL, PM, AL, AM, PC, C)

if(strcmpi( qd, 'pl'))
    x = PL(1);
    y = PL(2);
elseif (strcmpi( qd, 'pc'))
    x = PC(1);
    y = PC(2);
elseif (strcmpi( qd, 'pm'))
    x = PM(1);
    y = PM(2);
elseif (strcmpi( qd, 'am'))
    x = AM(1);
    y = AM(2);
elseif (strcmpi( qd, 'al'))
    x = AL(1);
    y = AL(2);
elseif (strcmpi( qd, 'c'))
    x = C(1);
    y = C(2);
else
    error('Wrong quadrant');
end;