function     [x,y]=get_orig(qd, PL, PM, AL, AM, PC, C)
% get the coordinates of an origin for a quadrant

if(strcmp( lower(qd), 'pl'))
    x = PL(1);
    y = PL(2);
elseif (strcmp( lower(qd), 'pc'))
    x = PC(1);
    y = PC(2);
elseif (strcmp( lower(qd), 'pm'))
    x = PM(1);
    y = PM(2);
elseif (strcmp( lower(qd), 'am'))
    x = AM(1);
    y = AM(2);
elseif (strcmp( lower(qd), 'al'))
    x = AL(1);
    y = AL(2);
elseif (strcmp( lower(qd), 'c'))
    x = C(1);
    y = C(2);
else
    disp('Wrong quadrant');
end;
