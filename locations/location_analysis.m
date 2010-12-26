
% in this file i wanted to find out, wether the neurons that are
% responsible for certain activation patterns have some topographical
% organization



conf.res_folder  = '~/Documents/uni/yifat_lab/results/';
conf.directory   = '/Volumes/LAB/';
conf.background  = 'mapgrid1bw.bmp';

% 1 for pronation supination, 2 for visual observable responses
conf.coloring           = 2;
conf.plot_connections   = 1;
conf.pro_sup_color      = {'r','b'};
conf.plot_field_size    = 1;
conf.monk               = 'vega';


% load a result file from the syn_analysis
load([conf.res_folder data filesep 'evoked_data_' conf.monk]);

% and load data about stimulation (origins of different reference frames)
load centers
load scale

% read the image in which to plot, paint it at set a origin of the plot with respect to the image.
curmap = imread(conf.background);
imagesc(curmap);
h = gcf;
colormap gray
set(h,'Position',[  1          29        1280         929]);
axis equal
hold on


% plot the origins in the map, with respect to which the stimulation sites were measured
plot(AL(1),AL(2),'rx');
plot(AM(1),AM(2),'rx');
plot(PL(1),PL(2),'rx');
plot(PM(1),PM(2),'rx');
plot(PC(1),PC(2),'rx');
plot(C(1), C(2), 'rx');


% add some information to the struct
for i = 1:length(resps)
    resps(i).plotted = 0; %#ok<SAGROW>
end

% for all sessions
for i=1:length(resps),
    
    % check whether the info file exists
    fname = [conf.directory conf.monk filesep 'info_files' filesep ...
        resps(i).sessions '_param.mat'];
    if ~exist(fname,'file'),
        disp(['Cannot find file --> ' fname]);
        
        % if info file for a session exists, load it and get stimulation
        % data.
    else
        load(fname);
                
        % TODO in vega I want to use only stimulation from electrode 2 and 3, 1 is spinal

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
                    if ~isempty(ctx_coord(1).qd) && ~isempty(ctx_coord(2).qd),
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
                
                sessdir(i).plotted = 1;
                sessdir(i).xy = [x y];
                
                if conf.coloring == 1
                    h = plot(x,y, ['o' conf.pro_sup_color{sessdir(i).hand}]);
                    set(h,'MarkerFaceColor',  conf.pro_sup_color{sessdir(i).hand});
                else
                    h = plot(x,y, ['o' ctx_coord(indx).eff]);
                    set(h,'MarkerFaceColor',  ctx_coord(indx).eff);
                end
                
                if conf.plot_field_size == 1
                    set(h,'MarkerSize', sessdir(i).field*2);
                end
                
            end
            
            % plot dashed connections between sites of two electrode stimulation
            if conf.plot_connections == 1
                if totel == 2,
                    plot(xx,yy,'k:');
                end
            end
        end
    end;
end
hold off


%%
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
end



figure
plot(loc.xy_dist, loc.ac_dist,'.');
[r p] = corrcoef(loc.xy_dist, loc.ac_dist)
figure
hist(loc.xy_dist)

figure
plot(loc.xy_dist(loc.xy_dist < 200), loc.ac_dist(loc.xy_dist < 200),'.');
[r p] = corrcoef(loc.xy_dist(loc.xy_dist < 200), loc.ac_dist(loc.xy_dist < 200))

