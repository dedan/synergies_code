function show_corticalmapchalva(path, monk)

srcdir = [path monk filesep];
load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);
load(['scales_' monk]);


% prepare the plot
figure(2)
clf
curmap = imread([monk '_cortex.bmp']);
imagesc(curmap);
colormap gray
axis equal
h = gcf;
set(h,'Position',[1 29 1280 929]);
hold on

if strcmp('vega', monk)
    plot(AL(1),AL(2),'rx');
    plot(AM(1),AM(2),'rx');
    plot(PL(1),PL(2),'rx');
    plot(PM(1),PM(2),'rx');
    plot(PC(1),PC(2),'rx');
    plot(C(1),C(2),'rx');
end


for i=1:length(resps),
    
    curdir   = char(resps(i).session);
    fname    = [curdir '_param.mat'];
    fullname = [srcdir 'info_files' filesep fname];
    
    if ~exist(fullname,'file'),
        disp(['Cannot find file --> ' fname]);
    else
        load(fullname);
        
        coord = new_cortical_data(DDFparam, SESSparam.SubSess(resps(i).subsession), monk);
        
        if strcmp(monk, 'chalva')
            
            % the 1.5 corresponde to the distance between the figure origin
            % to the actual in-chamber origin
            x = X0 - coord.x * D10/10;
            y = Y0 + coord.y * D10/10;
        elseif strcmp(monk, 'vega')
            % coordinate transformation of the data with respect to origin read from the
            % info file.
            if coord.posi == 0,         % means circular positioner
                [x,y]   = get_orig(coord.qd, PL, PM, AL, AM, PC, C);
                x       = x - coord.x / 5 * Scale5;
                y       = y - coord.y / 5 * Scale5;
                
            elseif coord.posi == 1,     % means dual positioner
                [x,y]   = get_orig('C', PL, PM, AL, AM, PC, C);
                [dx,dy] = translate_Dual(coord.electrode -1, curdir, coord.x, coord.y);
                x       = x + dx / 5 * Scale5;
                y       = y - dy / 5 * Scale5;
            else 
                error('bla');
            end
        end
        
        h = plot(x, y, 'ob');
                
        set(h,'MarkerSize', eps + resps(i).field *2);
        set(h,'MarkerFaceColor',  'b');
        
    end;
end
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord = new_cortical_data(ddf, subs, monk)

% which electrodes were used
% in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal
used_electrodes = find([ddf.Electrode.InUse]);
if strcmp(monk, 'vega') || strcmp(monk, 'vega_first')
    used_electrodes = used_electrodes(used_electrodes ~= 1);
end

for used = used_electrodes
    if isfield(subs.Electrode(used).Stim, 'Flag') && subs.Electrode(used).Stim.Flag
        coord.x         = ddf.Electrode(used).Y;
        coord.y         = ddf.Electrode(used).X;
        coord.posi      = ddf.Positioner;
        coord.qd        = ddf.Electrode(used).Quad;
        coord.electrode = used;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some translation with respect to date of recording session!??
function [dx,dy] = translate_Dual(electrode, curdir, x,y)

DD = str2double(curdir(2:3));
MM = str2double(curdir(4:5));

if (DD >= 15 && MM ==2) || (MM > 2),
    if electrode == 1, % left electrode
        dx = x + 6.5;
        dy = y + 1;
    elseif electrode == 2,
        dx = x - 4;
        dy = y - 2;
    end
elseif DD==14 && MM==2,
    if electrode == 1, % left electrode
        dx = x + 7.5;
        dy = y + .5;
    elseif electrode == 2,
        dx = x - 5;
        dy = y - 1;
    end
elseif DD>= 7 && MM== 2,
    if electrode == 1, % left electrode
        dx = x + 10.8;
        dy = y - 2;
    elseif electrode == 2,
        dx = x - 11.5;
        dy = y + 0.5;
    end
else
    error('dual translation problem');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
