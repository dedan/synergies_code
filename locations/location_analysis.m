
% in this file i wanted to find out, wether the neurons that are
% responsible for certain activation patterns have some topographical
% organization



config.res_folder  = '~/Documents/uni/yifat_lab/results/';
config.directory   = '/Volumes/LAB/vega/data/';
config.background  = 'mapgrid1bw.bmp';
config.plot_connections = 0;
config.plot_with_noise = 0;
config.coloring     = 1;        % 1 for pronation supination, 2 for visual observable responses
config.pro_sup_color = {'r','b'};
config.plot_field_size = 1;



% load a result file from the syn_analysis
load('latest_evoked');     %latest run of final code (I hope)

% and load data about stimulation (origins of different reference frames)
load centers
load scale

% read the image in which to plot, paint it at set a origin of the plot with respect to the image.
curmap = imread(config.background);
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
plot(C(1),C(2),'rx');


% put information which is available for the analysis into a struct
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
    indir = [config.directory curdir filesep];
    
    % check whether the info file exists
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
            
            % NOTE this is wrong. It was not stimulated with both
            % electrodes but the response is plotted on the location of
            % both electrodes.
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
                    
                    % put small noise on the coordinates that the dots don't cover eacht other
                    if config.plot_with_noise == 1
                        x = x + randn(1)*25;
                        y = y + randn(1)*25;
                    end
                    
                    if config.coloring == 1
                        cmnd = ['h=plot(x,y,''' 'o' config.pro_sup_color{sessdir(i).hand} ''');'];
                        eval(cmnd);
                        set(h,'MarkerFaceColor',  config.pro_sup_color{sessdir(i).hand});
                    else
                        cmnd = ['h=plot(x,y,''' 'o' ctx_coord(indx).eff ''');'];
                        eval(cmnd);
                        set(h,'MarkerFaceColor',  ctx_coord(indx).eff);
                    end
                    
                    if config.plot_field_size == 1
                        set(h,'MarkerSize', sessdir(i).field*2);
                    end
                    
                end
                
                % plot dashed connections between sites of two electrode stimulation
                if config.plot_connections == 1
                    if totel == 2,
                        plot(xx,yy,'k:');
                    end
                end
            end
        end;
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

