function show_mapvega( fin_res )

load centers
load scale

% load a result file from the syn_analysis
%load('~/Documents/uni/yifat_lab/results/data/120909_2012.mat');     %latest run of final code (I hope)

% directory_yifat = '\vega\data\';
directory = '/Volumes/DEDAN_DISK/vega/data/';


% read the image in which to plot, paint it at set a origin of the plot
% with respect to the image.
curmap = imread('mapgrid1bw.bmp'); 
imagesc(curmap);
h = gcf;
colormap gray
set(h,'Position',[  1          29        1280         929]);
axis equal
PenIndx =1;
hold on


% plot the origins in the map, with respect to which the stimulation sites
% were measured
plot(AL(1),AL(2),'rx');
plot(AM(1),AM(2),'rx');
plot(PL(1),PL(2),'rx');
plot(PM(1),PM(2),'rx');
plot(PC(1),PC(2),'rx');
plot(C(1),C(2),'rx');


% TODO output of a list of sessions used for this plot
slc = 0;

% get a list of all session folders
sessdir = dir([directory 'v*']);
syms = {'o','s'};


sessdir = struct;
for i = 1:length(fin_res(3).resp)
    sessdir(i).name = fin_res(3).resp(i).info.session;
    sessdir(i).field = fin_res(3).resp(i).field;
    sessdir(i).hand = fin_res(3).resp(i).hand;
    sessdir(i).resp = fin_res(3).dat(i,:);
    sessdir(i).plotted = 0;
end



% for all sessions
for i=1:length(sessdir),
    curdir = char(sessdir(i).name);
    indir = [directory curdir filesep];
    
    if exist([indir 'info'], 'dir'),
        fname = [curdir '_param.mat'];
        fullname = [indir 'info' filesep fname];
        if ~exist(fullname,'file'),
            disp(['Cannot find file --> ' fname]);
            
        % if info file for a session exists, load it and get stimulation
        % data.
        else
            load(fullname);
            ctx_coord = read_cortical_data( DDFparam);
            totel = 0;
            for indx = 1:2,
                if ~isempty(ctx_coord(indx).x),
                    totel = totel+1;
                    
                    % if there is not positioner flag in the info file, but there is 
                    % information about the quadrant for both electrodes availabe, positioner is 
                    % set to 1 (dual positioner);
                    if ~isfield(DDFparam, 'Positioner'),
                        disp(['No positioner flag in directory - ' curdir]);
                        if ~isempty(ctx_coord(1).qd) & ~isempty(ctx_coord(2).qd),
                             DDFparam.Positioner =1;
                        else
                             DDFparam.Positioner =0;
                        end;
                    end
                    
                    % coordinate transformation of the data with respect to origin read from the
                    % info file.
                    if DDFparam.Positioner ==0, % means circular positioner
                        [x,y]=get_orig(ctx_coord(indx).qd, PL, PM, AL, AM, PC, C);
                        x =x-ctx_coord(indx).x/5*Scale5;
                        y =y-ctx_coord(indx).y/5*Scale5;
                    elseif DDFparam.Positioner ==1, % means dual positioner
                        [x,y]=get_orig('C', PL, PM, AL, AM, PC, C);
                        [dx,dy] = translate_Dual( indx, curdir, ctx_coord(indx).x, ctx_coord(indx).y);
                        x =x+dx/5*Scale5;
                        y =y-dy/5*Scale5;
                        xx(indx)=x;
                        yy(indx)=y;
                    else 
                        disp('Stop!!!');
                    end;
                    
                    % chose a symbol which should be plotted with respect to the threshold value
                    if ctx_coord(indx).thr < 15,
                        symindx = 1;
                    else
                        symindx =2;
                    end
                    
                    sessdir(i).plotted = 1;
                    sessdir(i).xy = [x y];
                    
                    % put noise on the finally computed coordinates?
%                    x = x + randn(1)*25;
%                    y = y + randn(1)*25;


                    % plot the threshold chosen symbol at the coordinate.
                    
                    symindx = 1;
                    pro_sup_color = {'r','b'};
                    
                    cmnd = ['h=plot(x,y,''' syms{symindx} pro_sup_color{sessdir(i).hand} ''');'];
%                   cmnd = ['h=plot(x,y,''' syms{symindx} ctx_coord(indx).eff ''');'];
                     eval(cmnd);
                    
                    
                    
%                    set(h,'MarkerFaceColor',  ctx_coord(indx).eff);
                    set(h,'MarkerFaceColor',  pro_sup_color{sessdir(i).hand});
 %                   set(h,'MarkerSize', sessdir(i).field*2);
                    
                    cimenu(PenIndx) = uicontextmenu;
                    set(h,'UIContextMenu',cimenu(PenIndx) );
                    item1 = uimenu(cimenu(PenIndx), 'Label',[curdir '-' num2str(ctx_coord(1).id)]);
                end
                
                % plot dashed connections between sites of two electrode stimulation
                if totel == 2,
                    % plot(xx,yy,'k:');             
                end
            end
        end;
    end;
end
hold off


c = 0;
loc = struct;
for i = 1:length(sessdir)
    if sessdir(i).plotted == 1
        c = c+1;
        loc.xy(c,:)     = sessdir(i).xy;
        loc.resp(c,:)   = sessdir(i).resp;
    end
end


for i = 1:length(loc);
    loc.xy_dist      = pdist(loc.xy);
    loc.ac_dist      = pdist(loc.resp);
%    loc(i).range        = find(loc(i).xy_dist < 120);
end



figure
plot(loc.xy_dist, loc.ac_dist);
%plot_location(loc, 'ac_dist');








% get information about the cortical stimulation for a session
% information contains which quadrant, session id, number of electrodes
% used for stimulation, threshold and the response to this stimulation.
function ctx_coord = read_cortical_data( PRM)

    ctx_coord(1).id = PRM.ID;

for i=1:2,
    ctx_coord(i).qd = [];
    ctx_coord(i).x = [];
    ctx_coord(i).y = [];
    ctx_coord(i).x = [];
    ctx_coord(i).thr = [];
    ctx_coord(i).eff = [];
 end;
if isfield(  PRM, 'Cortex1'), % there are two electrodes
    if PRM.Cortex1.Flag
        ctx_coord(1).qd = PRM.Cortex1.Quad;
        ctx_coord(1).x = PRM.Cortex1.Y;
        ctx_coord(1).y = PRM.Cortex1.X;
        ctx_coord(1).thr = PRM.CTXmap(1).Thr;
        ctx_coord(1).eff = quantify_effect( PRM.CTXmap(1).Resp);
    end
    if PRM.Cortex2.Flag,
        ctx_coord(2).qd = PRM.Cortex2.Quad;
        ctx_coord(2).x = PRM.Cortex2.Y;
        ctx_coord(2).y = PRM.Cortex2.X;
        ctx_coord(2).thr = PRM.CTXmap(2).Thr;
        ctx_coord(2).eff = quantify_effect( PRM.CTXmap(2).Resp);
        
    else
        return;
    end;
else
    if PRM.Cortex.Flag,
        disp('Single electrode');
        ctx_coord(1).qd = PRM.Cortex.Quad;
        ctx_coord(1).x = PRM.Cortex.Y;
        ctx_coord(1).y = PRM.Cortex.X;
        ctx_coord(1).thr = PRM.CTXmap(1).Thr;
        ctx_coord(1).eff = quantify_effect( PRM.CTXmap(1).Resp);
    else
        return;
    end
end
            

% assign colors to responses of different body sites
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


% get the coordinates of an origin for a quadrant
function     [x,y]=get_orig(qd, PL, PM, AL, AM, PC, C)

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


% some translation with respect to date of recording session!??
function [dx,dy] = translate_Dual( iel, curdir, x,y);

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
    

        
        
