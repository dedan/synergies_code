function flat = flatten(grouped)

% i use this function with the group function 
% return of group is a struct and flatten can be used for example plotting
% the return with imagesc

flat = [];
for i = 1:length(grouped)
    flat = [flat; grouped(i).dat]; %#ok<AGROW>
end