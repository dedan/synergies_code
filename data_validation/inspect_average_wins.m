
%% debug window plots
% plot all the average windows, this can be used with the
% ../data_validation/sort_files.

% type == 1 is good to inspect whether an artefact reached into the window
% type == 2 to inspect the individual channels

function inspect_average_wins(path, monk, format, type)

folder_path = [path filesep 'data_validation' filesep 'average_wins'];


load([path 'data' filesep 'evoked_data_' monk])

% create subfolder if not exists or delete old files
if ~exist(folder_path, 'file')
    mkdir(folder_path);
else
    delete([folder_path filesep '*.' format]);
end


for i = 1:length(resps)
    
    % plot response 
    h = figure('Visible','off');
    tmp = resps(i).windows;
    
    if type == 1
        plot(resps(i).x, tmp');
        hold on
        plot(6*ones(1,11), 0:10, '--');
        legend(strtrim(cellstr(int2str(find(resps(i).c2take).'))));
        hold off
    else
        for j = 1:size(tmp,1)
            subplot(4,4,j);
            plot(resps(i).x, tmp(j,:));
        end
    end
    
    % indicate whether it is a session that was stored in several files
    if resps(i).connected
        title('connected');
    end
    
    saveas(h, [folder_path filesep monk(1) sprintf('%03d', i) '.' format]);
    close(h);
end

