function response_map(path, monk)

srcdir = [path monk filesep];
load([path 'results' filesep 'data' filesep 'evoked_data_' monk '.mat']);
load([path 'results' filesep 'data' filesep 'nat_mov_res_' monk '.mat']);
load(['scales_' monk]);


% prepare the plot
figure(2)
clf
curmap = imread([monk '_cortex.bmp']);
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


for i=1:length(resps),
    
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
                error('positioner variable not available');
            end
        end
        
        h = plot(x, y, 'ko');
        set(h,'MarkerSize', eps + resps(i).field *2);
        set(h,'MarkerFaceColor',  'b');
        
        c2take = nat_mov_res.c2take;
        [x1 y1] = pol2cart(nat_mov_res.pds(1,c2take), 1 *(resps(i).response - min(resps(i).response)));
        f = 1;
        plot([x, x + f *sum(x1)], [y, y + f * sum(y1)], 'k');
        
        
    end;
end
hold off


