function show_corticalmapchalva()

load('/Volumes/LAB/results/data/evoked_data_chalva.mat');
srcdir = '/Volumes/LAB/chalva/';
monk   = 'chalva';

figure(1)
clf

curmap = imread('chalva_cortex.bmp');
imagesc(curmap);
colormap gray
axis equal
h = gcf;
set(h,'Position',[  1          29        1280         929]);
load scales
hold on

dy = 1*(D10/10);
dx = 1*(D10/10);
plot(X0,Y0,'go');

% plotting a grid
for i=-10:10,
    for j=-10:10,
        cury = Y0 + dy*i;
        curx = X0 + j*dx;
        hh = plot(curx, cury, '.k');
        set(hh,'MarkerSize',1);
        set(hh,'Color',[.5 .5 .5]);
    end
end

sessdir = dir([srcdir filesep 'data' filesep monk(1) '*']);
syms    = {'o','^','s'};

PenIndx = 1;

for i=1:length(sessdir),
    
    curdir   = char(sessdir(i).name);
    fname    = [curdir '_param.mat'];
    fullname = [srcdir 'info_files' filesep fname];
    
    if ~exist(fullname,'file'),
        disp(['Cannot find file --> ' fname]);
    else
        load(fullname);
        sp_coord = read_cortical_data( DDFparam);
        
        
        % 
%        coord = new_cortical_data(DDFparam, sessions(i).subsession);
        
        
        curID    = DDFparam.ID;
        
        for j = 1:length(sp_coord)
            
            % the 1.5 corresponde to the distance between the figure origin
            % to the actual in-chamber origin
            x = X0-(sp_coord(j).x)*D10/10;
            y = Y0+(sp_coord(j).y)*D10/10;

            h = plot(x, y, 'ok');
            
            set(h,'MarkerFaceColor',  pass_color(sp_coord(j).pass));

            if sp_coord(j).Thr < 10 && sp_coord(j).Thr > 0,
                set(h,'MarkerSize',15);
            else
                set(h,'MarkerSize',10);
            end
            
            % context menu
            cimenu(PenIndx) = uicontextmenu;
            set(h,'UIContextMenu',cimenu(PenIndx) );
            uimenu(cimenu(PenIndx), 'Label',[curdir '-' num2str(curID)]);
            PenIndx = PenIndx + 1;
        end
    end;
end
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = pass_color(pass)
switch (pass),
    case 0,
        c = 'y';
    case 1,
        c = 'r';
    case 2,
        c = 'b';
    case 3,
        c = 'g';
    case 4,
        c = 'c';
    case 99,
        c = [0.8 0.8 0.8];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord = new_cortical_data(ddf, subs)

% which electrodes were used
% in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal
used_electrodes = find([DDFparam.Electrode.InUse]);
if strcmp(monk, 'vega') || strcmp(monk, 'vega_first')
    used_electrodes = used_electrodes(used_electrodes ~= 1);
end

for used = used_electrodes
    if isfield(subs.Electrode(used).Stim, 'Flag') && subs.Electrode(used).Stim.Flag
        coord.x     = DDFparam.Electrode(used).X;
        coord.y     = DDFparam.Electrode(used).Y;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctx_coord = read_cortical_data( PRM)

ctx_coord.x = [];
ctx_coord.y = [];
ctx_coord.id = PRM.ID;
indx = 1;
for i=1:length(PRM.Electrode),
    if PRM.Electrode(i).InUse,
        ctx_coord(indx).y = PRM.Electrode(i).Y;
        ctx_coord(indx).x = PRM.Electrode(i).X;
        ctx_coord(indx).pass = get_pass_type(PRM.Electrode(i).Active);
        ctx_coord(indx).Thr =  str2double(PRM.Electrode(i).Threshold);
        indx = indx+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   eff = quantify_effect( resp)

if ~isempty(findstr( lower(resp), 'finger'))
    eff = 'r';
elseif ~isempty(findstr( lower(resp), 'thumb'))
    eff = 'r';
elseif ~isempty(findstr( lower(resp), 'wrist'))
    eff = 'm';
elseif ~isempty(findstr( lower(resp), 'elbow'))
    eff = 'g';
elseif ~isempty(findstr( lower(resp), 'shoulder'))
    eff = 'b';
elseif ~isempty(findstr( lower(resp), 'sholder'))
    eff = 'b';
else
    disp(['Effect is --> ' resp]);
    eff = 'y';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function passtype = get_pass_type(PassResp)

passtype = 0;

disp(PassResp);
if isempty(PassResp),
    return;
end
PassResp = lower(PassResp);
if ~isempty(findstr(PassResp,'wrist')),
    passtype = 1;
elseif ~isempty(findstr(PassResp,'fingers')) | ~isempty(findstr(PassResp,'thumb')),
    passtype = 2;
elseif ~isempty(findstr(PassResp,'elbow')),
    passtype = 3;
elseif   ~isempty(findstr(PassResp,'shoulder')),
    passtype = 4;
else
    disp([PassResp '--> unidentified type'])
end


