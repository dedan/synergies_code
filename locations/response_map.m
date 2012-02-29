
data_path = '/Volumes/LAB/';
path = '~/projects/yifat_paper/';
monk = 'chalva';

srcdir = [data_path monk filesep];
load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);
load([path 'results' filesep 'data' filesep 'nat_mov_res_' monk '.mat']);
load(['data' filesep 'scales_' monk]);

dat = normc(vertcat(resps.response));
cluster_idx = kmeans(dat,4, 'replicates', 100);

colors = {'r', 'b', 'k', 'g'};
markers = {'^', 'v'};

% 1: no
% 2: cluster
% 3: 0-response
% 4: visible effect
coloring        = 3;
draw_direction  = false;
response_field  = false;
add_noise       = true;

% prepare the plot
f = figure('visible', 'off')
clf
curmap = imread(['data' filesep monk '_cortex.bmp']);
imagesc(curmap);
colormap gray
axis     equal
hold     on
axis([0 size(curmap,2) 0 size(curmap,1)])

if strcmp('vega', monk)
    plot(AL(1),AL(2),'rx');
    plot(AM(1),AM(2),'rx');
    plot(PL(1),PL(2),'rx');
    plot(PM(1),PM(2),'rx');
    plot(PC(1),PC(2),'rx');
    plot(C(1),C(2),'rx');
end


for i=1:length(resps)

    curdir   = char(resps(i).session);
    fname    = [curdir '_param.mat'];
    fullname = [srcdir 'info_files' filesep fname];

    if ~exist(fullname,'file'),
        disp(['Cannot find file --> ' fname]);
    else
        load(fullname);

        coord = get_cortical_data(DDFparam, SESSparam.SubSess(resps(i).subsession), monk);

        if strcmp(monk, 'chalva')

            % the 1.5 corresponde to the distance between the figure origin
            % to the actual in-chamber origin
            x = X0 - coord.x * D10/10;
            y = Y0 + coord.y * D10/10;
            plot(X0,Y0,'go');

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
                y       = y + dy / 5 * Scale5;
            else
                error('positioner variable not available');
            end
        end

        if add_noise
            x = x + randn(1)*4;
            y = y + randn(1)*4;
        end

        h = plot(x, y, 'ko');

        if response_field
            set(h,'MarkerSize', eps + resps(i).field *2);
        end

        if coloring == 1
            set(h,'MarkerFaceColor',  'none');
        elseif coloring == 2
            set(h,'MarkerFaceColor',  colors{cluster_idx(i)});
        elseif coloring == 3
            idx = (resps(i).field == 0) +1;
            set(h,'Marker',  markers{idx});
            set(h,'MarkerEdgeColor',  colors{idx});
            set(h,'MarkerSize', 7);
        elseif coloring == 4
            set(h,'MarkerEdgeColor',  quantify_effect(resps(i).location.res));
        end

        if draw_direction
            c2take = nat_mov_res.c2take;
            [x1 y1] = pol2cart(nat_mov_res.pds(1,c2take), 1 *(resps(i).response - min(resps(i).response)));
            f = 1;
            plot([x, x + f *sum(x1)], [y, y + f * sum(y1)], 'k');
        end


    end;
end
hold off

if coloring == 4
    disp('fingers -> r');
    disp('wrist -> g');
    disp('elbow -> b');
    disp('shoulder -> m');
    disp('face, body, none -> k');
end

saveas(f, [path 'results' filesep 'maps' filesep 'map_' monk '.png']);


