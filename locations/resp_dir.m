
% In this file I started search for something like an "orientation map" on
% the cortex. The Idea was to look at the Torque which was induced by the
% stimulation. When this is found to be consistent, not uniform or
% something else interesting, we could plot a map of the cortex with
% preferred directions of stimulation sites.

% no usefull results from this analysis yet

addpath(genpath('..'));
config.stim_value          = 250;                                          % look only at stimulations around this value
config.monkey              = 'vega';                                       % monkey to analyze
config.volume              = '/Volumes/lab/';                              % harddisk which contains the data
config.folder              = [config.volume config.monkey filesep];
config.pd_folder           = [config.folder 'pd_files' filesep];           % PD files (preferred direction)
config.dat_folder          = [config.folder 'data' filesep];               % folder with data


% get all stimulations
data = get_all_stimulations(config);

% filter subessions according to StimAmp
% here I look at stronger stimulations as in the syn_analysis as I hoped to
% get stronger Torque response from it
data = stimulations_at(data, config.stim_value);

i = 1;
valid_data = struct;
for dati = 1:length(data)
    
    dat = data(dati);
    for file_number = dat.file(1):dat.file(end)
        file = [config.dat_folder dat.session '/MAT/' dat.session lead_zeros(file_number,3) int2str(file_number)];
        
        try
            % check wether a StimTime vector exists
            stim_check = whos('-file',[file '_bhv'],'StimTime');
            
            % this filters out sessions which only have a few stimulations
            % (occurs when a file was divided at end of recording)
            if(stim_check.size(1) > 4)
                valid_data(i).file = file;
                i = i+1;
            end
        catch err
            display(['no StimTime Entry for: ' file ]);
        end
    end
end


%%
c = 1;
res = struct;
for i = 1:length(valid_data)
    %for i = 1
    file_base = valid_data(i).file;
    load([file_base '_bhv'], 'TrqFE', 'TrqRU', 'StimTime');
    try
        for j = 1:length(StimTime)
            %        plot(TrqFE(floor(StimTime(j)*1000 - 100):floor(StimTime(j)*1000 + 100)));
            if(floor(StimTime(j)*1000 + 100) < length(TrqFE) && floor(StimTime(j)*1000 - 100) > 0);
                range = floor(StimTime(j)*1000 - 100):floor(StimTime(j)*1000 + 100);
                if size(range,2) > 201
                    range = range(1:201);
                    c = c+1;
                end
                res(i).avgf(j,:) = TrqFE(range);
                res(i).avgf(j,:) = res(i).avgf(j,:) - mean(res(i).avgf(j,1:100));
                res(i).avgu(j,:) = TrqRU(range);
                res(i).avgu(j,:) = res(i).avgu(j,:) - mean(res(i).avgu(j,1:100));
                res(i).x(j)     = sum(res(i).avgf(j,106:120)) / sqrt(sum(res(i).avgf(j,106:120))^2 + sum(res(i).avgu(j,106:120))^2);
                res(i).y(j)     = sum(res(i).avgu(j,106:120)) / sqrt(sum(res(i).avgf(j,106:120))^2 + sum(res(i).avgu(j,106:120))^2);
            end
        end
    catch err
        disp(err)
    end
end

c

%%
for i = 1:length(res)
    subplot 311
    plot(-100:100,mean(res(i).avgf));
    subplot 312
    plot(-100:100,mean(res(i).avgu));
    subplot 313
    compass(res(i).x, res(i).y);
    pause;
end



