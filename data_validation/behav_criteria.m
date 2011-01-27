
% previously the selection of valid trials was done just as yuval did it in
% his study without questioning this criteria again. This script
% investigates the distribution of variables which are used as selection
% criterial. We decided not to be very strict in the selection but just to
% remove the long tail.

% STEPHAN: here we select trials based on some behavioral
% criteria of reaction time, {-200-500) Movement time
% (500-1500) and some angular deviation (<35). I think we
% could relax some of these criteria) especially the last
% one

% path = 'E:\';
clear
base_path    = '/Volumes/LAB/';
outpath = '~/Documents/uni/yifat_lab/results/data_validation/';
names   = {'vega', 'darma', 'chalva'};
e2add   = 'e';
fields  = {'reaction', 'movement', 'ang_div'};
data    = struct;

set.reaction_win = [-200, 500];
set.movement_win = [0, 1500];
set.ang_div_win  = [-35, 35];

set.reaction_lim = [-1000,2000];
set.movement_lim = [-1000,2000];
set.ang_div_lim  = [-100,100];


if exist([outpath 'criteria_dist.mat'], 'file')
    load([outpath 'criteria_dist']);
else
    for l = 1:length(names)
        
        data(l).reaction = [];
        data(l).movement = [];
        data(l).ang_div  = [];
        
        monk = char(names{l});
        path     = [base_path monk filesep];
        load([path monk 'session']);
        
        vdir = dir([path 'data' filesep monk(1) '*']);
        vdir = sortdirs( vdir);
        
        for i= 1:length(vdir),
            
            sessname = char(vdir(i));
            disp(['file ' num2str(i) ' of ' num2str(length(vdir))]);
            
            emgdir = [path 'data' filesep sessname filesep 'mat' filesep ];
            eddir  = [path 'EDfiles' filesep];
            
            femgs = dir([emgdir '*_emg.mat']);
            if isempty(femgs),
                disp('No EMG mat files');
                f1      = 0;
                f2      = 0;
            else
                f1 = 1;
                f2 = length(femgs);
            end
            
            
            indx        = 1;
            point2empty = [];
            edfiles     = cell(0);
            
            for j=f1:f2,
                edname = extract_edname( monk(1), subss, sessname,  j, e2add );
                if ~isempty(edname),
                    edfiles(j-f1+1)     = {[eddir edname]};
                else
                    point2empty(indx)   = j-f1+1;  %#ok<SAGROW>
                    indx                = indx+1;
                end
            end
            
            if isempty(edfiles),
                disp('empty edfiles');
            end
            
            if ~isempty(point2empty),
                for j=1:length(point2empty),
                    curpos = point2empty(j);
                    disp('--> no edfiles was found');
                end
                indx = ones(size(f1:f2));
                indx(point2empty) = 0;
                edfiles = edfiles(logical(indx));
            end
            
            for j=1:length(edfiles)
                
                if exist(char(edfiles(j)), 'file')
                    bhvdata = load(char(edfiles(j)));
                    if isfield(bhvdata, 'bhvStat') && ~isempty(bhvdata.bhvStat)
                        bhvStat = bhvdata.bhvStat;
                        
                        data(l).reaction = [data(l).reaction; bhvStat(:,1)];
                        data(l).movement = [data(l).movement; bhvStat(:,2)];
                        data(l).ang_div  = [data(l).ang_div;  bhvStat(:,5)]; %#ok<*AGROW>
                    end
                end
            end
        end
    end
    
    save([outpath 'criteria_dist'], 'data');
end


for i = 1:length(fields)
    
    h = figure('Visible', 'off');
    
    for j = 1:length(names)
        subplot(length(names)+1, 1, j)
        
        hist(data(j).(fields{i}), 100)
        hold on
        xlim(set.([fields{i} '_lim']));
        tmp = set.([fields{i} '_win']);
        x1  = ones(1,10) * tmp(1);
        x2  = ones(1,10) * tmp(2);
        y   = 1:500:5000;
        plot(x1, y, 'r.');
        plot(x2, y, 'r.');
        hold off
        title([fields{i} ' - ' names{j}]);
        
    end
    subplot(length(names)+1, 1, length(names)+1)
    
    hist(vertcat(data.(fields{i})), 100)
    hold on
    xlim(set.([fields{i} '_lim']));
    tmp = set.([fields{i} '_win']);
    x1  = ones(1,10) * tmp(1);
    x2  = ones(1,10) * tmp(2);
    y   = 1:500:5000;
    plot(x1, y, 'r.');
    plot(x2, y, 'r.');
    hold off
    title([fields{i} ' - all']);
    saveas(h, ['~/Documents/uni/yifat_lab/results/data_validation/' fields{i} '.pdf']);
    close(h);
end