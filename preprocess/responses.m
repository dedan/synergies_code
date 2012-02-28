function resp = responses(filtered_data, config)

% get all the responses (integrated average windows) for all sessions.
% as we have only very little session in the midposition posture, we ignore
% them at all

resp = struct;

j = 1;
for dati = 1:length(filtered_data)

    if 1 %(config.inflag == 1)
        disp(' ');
        display(['working on session ' int2str(dati) ' of ' int2str(length(filtered_data))]);
    end

    dat = filtered_data(dati);


    % collect the data for the average window function
    emg_data = [];
    StimTime_agg = [];

    % aggregate divided files
    if length(dat.file(1):dat.file(end)) > 1
        disp('have to connect');
    end

    files = cell(1,length(dat.file(1):dat.file(end)));
    for i = 1:length(dat.file(1):dat.file(end))
        dat.file = dat.file(1):dat.file(end); % YP there was a bug that happened rarely in cases were we had 3 stim files.
        file_number = dat.file(i);
        files{i}    = [dat.session sprintf('%03d',file_number)];
        file        = [config.dat_folder 'data' filesep dat.session ...
                        filesep 'MAT' filesep files{i}];

        stim_check = whos('-file',[file '_bhv'],'AMstim_on', 'StimTime');
        if(stim_check.size(1) > 4) && exist([file '_emg.mat'], 'file')


            s_times = load([file '_bhv.mat'], stim_check.name);
            s_times = s_times.(stim_check.name);
            emg = load([file '_emg.mat'], 'EMG*');
            first_good = find(config.c2take,1);
            f_orig = emg.(['EMG' int2str(config.channels(first_good)) '_KHz']);

            tmp = NaN(length(find(config.c2take)), ...
                length(emg.(['EMG' int2str(config.channels(first_good))])));
            c = 1;
            for i = find(config.c2take)
                tmp(c,:) = emg.(['EMG' int2str(config.channels(i))]);
                c = c+1;
            end

            s_times = s_times +(size(emg_data,2) / (f_orig*1000));
            StimTime_agg = [StimTime_agg; s_times]; %#ok<AGROW>
            emg_data = [emg_data tmp]; %#ok<AGROW>


        end
    end


    if length(StimTime_agg) > 1

        % get hand information from ed file
        ed_name = [config.monk(1) sprintf('%02d', dat.id) sprintf('%02d', dat.subsession) ...
                    'ee.1.mat'];
        hand    = load([config.dat_folder 'EDfiles' filesep ed_name], 'hand_position');
        resp(j).hand        = hand.hand_position;
        resp(j).f_orig      = f_orig;
        resp(j).id          = dat.id;
        resp(j).session     = dat.session;
        resp(j).subsession  = dat.subsession;
        resp(j).amp         = dat.amp;
        resp(j).electrode   = dat.electrode;
        resp(j).location    = dat.location;
        resp(j).files       = files;
        resp(j).connected   = length(dat.file(1):dat.file(end)) > 1;

        % NOTE continue to skip the non pronation or supination sessions
        if(resp(j).hand == 3)
            if config.inflag
                display('hand in midposition, skip this session..');
            end
            continue;
        end

        wins = get_average_windows(emg_data, StimTime_agg, f_orig, config);

        resp(j).x       = wins.x_axis;
        resp(j).windows = wins.windows;

        % determine fieldsize (fieldsize is the number of
        % significantly responding channels)
        resp(j).field   = length(find(wins.p < 0.05 & wins.p > 0));
        resp(j).fields  = find(wins.p < 0.05 & wins.p > 0);

        for i = 1:size(wins.windows,1)

            % response computation as described by yuval
            mean_background = mean(resp(j).windows(i,wins.pre_r));
            mean_response = mean(resp(j).windows(i,wins.post_r));
            resp(j).response(i) = mean_response / mean_background;
        end

        j = j+1;
    else
        disp('session skipped');
    end
end


