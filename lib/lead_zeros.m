function res = lead_zeros(num, len)

% can be use, if for example in a file name a certain space for the number
% is reserved. So lead_zeros adds leading numbers to an int in order to get
% the right size

res = [];
while length(int2str(num)) < len
    res = [res '0']; %#ok<AGROW>
    num = num*10;
end
