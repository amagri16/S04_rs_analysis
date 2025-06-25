%% TIAD2 EEG Preprocessing

% Last edited by S. Williamson, 16/12/24

% Current sections
  % - Load in Data
  % - Define the task
  % - Low-pass filter
  % - Extract markers and downsample
  % - DFT and High-pass Filter
  % - FFT
  % - Check number of trials

  function Preproc(in_direc, out_direc, subject, downsample_freq, check_fft, verb)
    
    %% Loop through all files:
    
    for week = 1:2
        disp('week')
        disp(week)
    
        for study_day = 1:8
            disp('study_day')
            disp(study_day)
    
            if study_day > 1 && study_day < 6
                sessions = 2;
            else
                sessions = 1;
            end
    
            for session = 1:sessions
                disp('session')
                disp(session)
    
                if study_day > 1 && study_day < 6 && (study_day == 4 && session == 1) == 0
                    % stim tasks
                    tasks = 10;
                else 
                    % pre / post tasks
                    tasks = 13;
                end
    
                for task_no = 1:tasks
                    disp('task_no')
                    disp(task_no)
    
                    %% Set Parameters

                    if study_day == 1
                        day = 'pre_day';
                    end
                    if study_day == 2
                        day = 'stim_day1';
                    end
                    if study_day == 3
                        day = 'stim_day2';
                    end
                    if study_day == 4
                        day = 'stim_day3';
                    end
                    if study_day == 5
                        day = 'stim_day4';
                    end
                    if study_day == 6
                        day = 'post_day';
                    end
                    if study_day == 7
                        day = 'follow_up1';
                    end
                    if study_day == 8
                        day = 'follow_up3';
                    end
    
                    if study_day > 1 && study_day < 6 && (study_day == 4 && session == 1) == 0
                        if task_no == 1
                            task = 'rs';
                            run_no = '1';
                        end
                        if task_no == 2
                            task = 'f';
                            run_no = '1';
                        end
                        if task_no == 3
                            task = 'f';
                            run_no = '2';
                        end
                        if task_no == 4
                            task = 'wwa';
                            run_no = '1';
                        end
                        if task_no == 5
                            task = 'wwa';
                            run_no = '2';
                        end
                        if task_no == 6
                            task = 'pal';
                            run_no = '1';
                        end
                        if task_no == 7
                            task = 'pal';
                            run_no = '2';
                        end
                        if task_no == 8
                            task = 'osa';
                            run_no = '1';
                        end
                        if task_no == 9
                            task = 'osa';
                            run_no = '2';
                        end
                        if task_no == 10
                            task = 'rs';
                            run_no = '2';
                        end
    
                    else
    
                        if task_no == 1
                            task = 'rs';
                            run_no = '1';
                        end
                        if task_no == 2
                            task = 'om';
                            run_no = '';
                        end
                        if task_no == 3
                            task = 'wm';
                            run_no = '';
                        end
                        if task_no == 4
                            task = 'fnat';
                            run_no = '1';
                        end
                        if task_no == 5
                            task = 'fnat';
                            run_no = '2';
                        end
                        if task_no == 6
                            task = 'fnat';
                            run_no = '3';
                        end
                        if task_no == 7
                            task = 'om';
                            run_no = 'd';
                        end
                        if task_no == 8
                            task = 'wm';
                            run_no = 'd';
                        end
                        if task_no == 9
                            task = 'mst';
                            run_no = '1';
                        end
                        if task_no == 10
                            task = 'mst';
                            run_no = '2';
                        end
                        if task_no == 11
                            task = 'mst';
                            run_no = '3';
                        end
                        if task_no == 12
                            task = 'fnat';
                            run_no = 'd';
                        end
                        if task_no == 13
                            task = 'rs';
                            run_no = '2';
                        end
                    end
                    
                   
                    
                    %% Confirm output directory:

                    % Create the full save directory path
                    save_dir = fullfile(out_direc, 'preprocessed', ['subject_', subject], task, ['week_', num2str(week)]);

                    % Check if the directory exists, and create it if not
                    if ~exist(save_dir, 'dir')
                        mkdir(save_dir);
                    end

                    % Create the filename within the save directory
                    f_name = sprintf('%s_Session_%s_%s%s.mat', day, num2str(session), task, run_no);
                    full_path = fullfile(save_dir, f_name);

                    % Check if the file exists, and skip loop if it does
                    if exist(full_path, 'file')
                        fprintf('File %s already exists. Skipping save.\n', full_path);
                        continue
                    end

                    
                    %% Load in data:
                              
                    %path = [in_direc, filesep, 'subject_', subject, filesep, 'EEG', ...           % Create path
                    %    filesep 'week_', num2str(week), filesep, day];

                    path = fullfile(in_direc, ['subject_', subject], 'EEG', ['week_', num2str(week)], day);
    
                    % Check if sessions are present, extend path if so
                    % (session folders don't exist on days where there is only one session)
                    if ~isempty(dir(fullfile(path, ['session', '*'])))
                        path = [path, filesep, 'session_', num2str(session)];
                    end 
    
                    file = [task, run_no, '.vhdr'];                                            % Define file
    
                    disp([path, ' --> ', file])
                    EEG = pop_loadbv(path, file);                                              % Load file
                    
                    
                    if matches(task, 'mst') || matches(task, 'fnat') || matches(task, 'f')     % Determine if cognitron task or not
                        disp("non-cognitron task confirmed")
                        cognitron = false;
                    elseif matches(task, 'wm') || matches(task, 'om') || matches(task, 'wwa') || matches(task, 'pal') || matches(task, 'osa')
                        disp("cognitron task confirmed")
                        cognitron = true;
                    elseif matches(task, 'rs')
                        disp("resting state recording confirmed")
                        cognitron = false;
                    else
                        error("error: task type not found")
                    end 
                    
                    n_chan = EEG.nbchan - 1;                                                   % Define number of EEG channels
                    
                    %% Low-pass filter the data:
                    
                    EEG.data = EEG.data/10;                                                    % Undo amplification 
                    EEG = pop_eegfiltnew(EEG, 'hicutoff',downsample_freq/4,'plotfreqz',0);     % Set filter cutoff to half the nyquist freq. of downsampled sampling rate
                    
                    %% Extract markers and Downsample:
                    
                    if cognitron == true
                    
                        EEG = pop_resample(EEG, downsample_freq);                              % Downsample first if cognitron task
                        LED_signal = EEG.data(n_chan+1, :);                                    % Isolate the LED signal
                        LED_deriv1 = [0, diff(LED_signal)];                                    % Find 1st derivative peaks and troughs
                        
                        figure;
                        plot(LED_deriv1)
                        title([path, ' --> ', file])
                                                                                               % No of peaks to be found depends on no of trials
                        task_peak_no = containers.Map;                                         % Define expected trial numbers
                        task_peak_no('om') = 20;                                               % applies to omd too
                        task_peak_no('wm') = 24;                                               % applies to wmd too
                        task_peak_no('wwa') = [15, 18];
                        task_peak_no('pal') = 1;
                        task_peak_no('osa') = 18;

                        temp = task_peak_no(task);
                        max_peaks_allowed = temp(end);
                        min_peaks_allowed = temp(1);

                        if matches(task, 'pal')

                            alpha = 1;
                            peaks = ones(1, 100);                                              % placeholder length
                            while length(peaks) > 50
    
                                threshold = mean(LED_deriv1) + alpha*std(LED_deriv1); 
                                [peaks, peak_times] = findpeaks(LED_deriv1, ...
                                    'MinPeakHeight', threshold);
    
                                alpha = alpha + 1;
                            end 
                            diffs = peak_times(2:end) - peak_times(1:end-1);
                            diffs = [51, diffs];
                            peaks = peaks(diffs>50);
                            peak_times = peak_times(diffs>50);
                            if sum(diffs<50) > 0
                                 peaks = peaks(1:end-1);    
                                 peak_times = peak_times(1:end-1);
                            end
                            task_peak_no('pal') = length(peaks);

                        else

                            alpha = 10;
                           
                            peak_heights = [];
                            while length(peak_heights) < min_peaks_allowed
    
                                threshold = mean(LED_deriv1) + alpha*std(LED_deriv1);          % Set a reasonable threshold for finding peaks
    
                                [peak_heights, peak_times] = findpeaks(LED_deriv1, ...
                                    'MinPeakHeight', threshold);
                                diffs = peak_times(2:end) - peak_times(1:end-1);               % Here we are finding peaks, excluding spurious ones, counting peaks
                                diffs = [51, diffs];                                           %, and if there are fewer than expected we are lowering the threhsold and 
                                if length(peak_heights) >= min_peaks_allowed
                                    spurious_height_threshold = mean(peak_heights(1:6))*2;
                                    spurious_time_threshold = 50;
                                    spurious_peaks = (diffs<spurious_time_threshold) + (peak_heights>spurious_height_threshold);
                                    if sum(spurious_peaks) > 0
                                        indices = find(spurious_peaks > 0);
                                        first_spurious = indices(1);
                                        spurious_lower_bound = peak_times(first_spurious) - EEG.srate*3;
                                    end
                                    if sum(spurious_peaks) > 0 % && first_spurious > min_peaks_allowed
                                        peak_heights = peak_heights(1:first_spurious-1);                                   % trying again until we find the expected number. 
                                        peak_times = peak_times(1:first_spurious-1); 
                                        while peak_times(end) > spurious_lower_bound
                                            peak_heights = peak_heights(1:end-1);                                   % trying again until we find the expected number. 
                                            peak_times = peak_times(1:end-1);
                                        end
                                    end
                                    
                                    if length(peak_heights) > max_peaks_allowed     
                                        difference = length(peak_heights) - task_peak_no(task);
                                        peak_heights = peak_heights(1:end-difference);
                                        peak_times = peak_times(1:end-difference);
                                    end
                                end
                               
                                alpha = alpha - 1;
    
                                if alpha == 1
                                    error("error: unable to find peaks in the LED signal")
                                end 
    
                            end 
                        end


                        troughs = [];
                        while length(troughs) < length(peak_heights)

                            threshold = mean(LED_deriv1) + alpha*std(LED_deriv1);      
                            
                            LED_deriv1_inverted = -LED_deriv1;                                     % plot(LED_deriv1)  % plot(LED_deriv1_inverted)
                            [troughs, trough_times] = findpeaks(LED_deriv1_inverted, ...
                                'MinPeakHeight', threshold);
                            diffs = trough_times(2:end) - trough_times(1:end-1);
                            diffs = [51, diffs];
                            troughs = troughs(diffs>50);
                            trough_times = trough_times(diffs>50);
                            if sum(diffs < 50) > 0
                                troughs = troughs(1:end-1);
                                trough_times = trough_times(1:end-1);
                            end
                            if length(troughs) > length(peak_heights)                                     % if troughs are more than peaks after removal of spurious ones, assume the final troughs are spurious
                                difference = length(troughs) - length(peak_heights);
                                troughs = troughs(1:end-difference);
                                trough_times = trough_times(1:end-difference);
                            end
                            if length(troughs) < length(peak_heights)
                                difference = length(peak_heights) - length(troughs);
                                if difference >=2                                   % we don't want unreasonable subtractions
                                    alpha = alpha -1;
                                    continue
                                end 
                                peak_heights = peak_heights(1:end-difference);
                                peak_times = peak_times(1:end-difference);
                            end
                            alpha = alpha -1;

                            if alpha == 1
                                error("error: unable to find troughs in the LED signal");
                            end 

                        end 
                      
                        figure;
                        plot(1, 1)
                        title(['number extracted = ', num2str(length(peak_heights))])
                        
                        disp("Adding optical events to EEG data structure")                    % Add markers to EEG structure 
                        for i = 1:length(trough_times)                                         % Add each trial in turn (peaks are responses, troughs are stimuli / trial starts)
                          disp(i);
                          numEvents = length(EEG.event);                                       % Add trial start event:
                          EEG.event(numEvents + 1).type = 'S  2';                              % 'S2' so as to match non-cognitron tasks
                          EEG.event(numEvents + 1).latency = trough_times(i);
                          numEvents = length(EEG.event);                                       % Add response event:
                          EEG.event(numEvents + 1).type = 'response';
                          EEG.event(numEvents + 1).latency = peak_times(i);
                        end
                    
                    end 
                    
                    
                    if cognitron == false 
                    
                        if matches(task, 'mst')
                    
                            disp('No marker extraction steps required');
                            if verb                                                            % Print markers, latencies, and successive time differences
                                events = [EEG.event.latency];                                      
                                types = [EEG.event.type];
                                type_split = strsplit(types, ' ');             
                                events = events/EEG.srate;     
                                combined = [num2cell(events)', type_split'];                   %#ok<NASGU> 
                                event_diff = [0, diff(events)];
                                combined = [num2cell(events)', type_split', ...
                                    num2cell(event_diff)'];
                                disp(combined);
                            end 
                        end
                    
                    
                        if matches(task, 'fnat')
                    
                            tmp_lab = {EEG.event.type};                                        % Convert EEG.event.type to cell_array
                            is_S1 = strcmp(tmp_lab, 'S  1');                                   % Create logical arrays to find matches
                            is_S2 = strcmp(tmp_lab, 'S  2');
                            is_S3 = strcmp(tmp_lab, 'S  3');
                            is_S4 = strcmp(tmp_lab, 'S  4');
                            is_S5 = strcmp(tmp_lab, 'S  5');
                            is_S6 = strcmp(tmp_lab, 'S  6');
                            is_S7 = strcmp(tmp_lab, 'S  7');                                   % Remove spurious markers:
                            to_remove_1 = [false, false, is_S1(3:end)];                        % When 'S  1' repeats
                            to_remove_2_2 = [false, is_S2(2:end) & is_S2(1:end-1)];            % When 'S  2' follows 'S  2'
                            to_remove_4_2 = [false, is_S2(2:end) & is_S4(1:end-1)];            % When 'S  2' follows 'S  4'
                            to_remove_5_4 = [false, is_S4(2:end) & is_S5(1:end-1)];            % When 'S  4' follows 'S  5'
                            to_remove_6_4 = [false, is_S4(2:end) & is_S6(1:end-1)];            % When 'S  4' follows 'S  6'
                            to_remove_7_4 = [false, is_S4(2:end) & is_S7(1:end-1)];            % When 'S  4' follows 'S  7'
                            to_remove_4_6 = [is_S4(1:end-1) & is_S6(2:end), false];            % When 'S  4' precedes 'S  6'
                            to_remove_2_6 = [is_S2(1:end-1) & is_S6(2:end), false];            % When 'S  2' precedes 'S  6'
                            to_remove_3_7 = [is_S3(1:end-1) & is_S7(2:end), false];            % When 'S  3' precedes 'S  7'
                            to_remove_2_4 = [is_S2(1:end-1) & is_S4(2:end), false];            % When 'S  2' precedes 'S  4'
                            to_remove = to_remove_1 | to_remove_2_2 | ...
                                        to_remove_4_2| to_remove_5_4 | ...
                                        to_remove_6_4 | to_remove_7_4 | ...
                                        to_remove_4_6 | to_remove_2_6 | ...
                                        to_remove_3_7 | to_remove_2_4;
                            EEG.event(to_remove) = [];                                         % Remove markers based on the logical indices 
                            disp([num2str(sum(to_remove)), ' spurious markers removed']);      % Indicate how many were removed
                        
                            events = [EEG.event.latency];                                      % Complete latency-based marker removal also          
                            events = 1000*(events/EEG.srate);                                  % Convert latencies into milliseconds
                            threshold = 20;                                                    % Set threshold for removal of markers e.g. 10ms 
                            event_diff = [threshold, diff(events)];                            % Find successive time differences between markers
                            reasonably_delayed_markers = event_diff >= threshold;              % Remove any markers that follow previous markers by < threshold (ms)
                            markers_to_remove = event_diff < threshold;  
                            tmp_lab = {EEG.event.type};   
                            is_S2 = strcmp(tmp_lab, 'S  2');                                   % Recalculate S2 and S4 markers so that none are removed
                            is_S4 = strcmp(tmp_lab, 'S  4');
                            markers_to_remove(is_S2) = 0;
                            markers_to_remove(is_S4) = 0;
                            EEG.event(markers_to_remove) = [];
                            disp([num2str(sum(markers_to_remove)), ...                         % Indicate how many markers removed
                               ' markers removed due to following the previous marker by <' ...
                               , num2str(threshold), 'ms'])
                            
                            
                            if verb                                                            % Print remaining markers, latencies, and successive time differences
                                events = [EEG.event.latency];         
                                types = [EEG.event.type];
                                type_split = strsplit(types, ' '); 
                                events = events/1000;
                                combined = [num2cell(events)', type_split'];                   %#ok<NASGU> 
                                event_diff = [0, diff(events)];
                                combined = [num2cell(events)', type_split', ...
                                    num2cell(event_diff)'];
                                disp(combined);
                            end 
                        end 
                    
                    
                        if matches(task, 'f')
                    
                            tmp_lab = {EEG.event.type};
                            is_S1 = strcmp(tmp_lab, 'S  1');
                            is_S2 = strcmp(tmp_lab, 'S  2');
                            is_S3 = strcmp(tmp_lab, 'S  3');
                            is_S4 = strcmp(tmp_lab, 'S  4');
                            is_S5 = strcmp(tmp_lab, 'S  5');
                            is_S6 = strcmp(tmp_lab, 'S  6');
                            is_S7 = strcmp(tmp_lab, 'S  7');
                            is_S8 = strcmp(tmp_lab, 'S  8');
                            is_S10 = strcmp(tmp_lab, 'S 10');
                            is_8R = strcmp(tmp_lab, 'R  8');
                            to_remove_1 = [false, false, is_S1(3:end)];                        % When 'S  1' repeats            
                            to_remove_2_3 = [false, is_S2(2:end) & is_S3(1:end-1)];            % When 'S  2' follows 'S  3'      
                            to_remove_3_4 = [false, is_S3(2:end) & is_S4(1:end-1)];            % When 'S  3' follows 'S  4'    
                            to_remove_5_4 = [false, is_S4(2:end) & is_S5(1:end-1)];            % When 'S  4' follows 'S  5'        
                            to_remove_6_4 = [false, is_S4(2:end) & is_S6(1:end-1)];            % When 'S  4' follows 'S  6'        
                            to_remove_7_4 = [false, is_S4(2:end) & is_S7(1:end-1)];            % When 'S  4' follows 'S  7'       
                            to_remove_7_6 = [false, is_S6(2:end) & is_S7(1:end-1)];            % When 'S  6' follows 'S  7'      
                            to_remove_2_6 = [is_S2(1:end-1) & is_S6(2:end), false];            % When 'S  2' precedes 'S  6'     
                            to_remove_3_7 = [is_S3(1:end-1) & is_S7(2:end), false];            % When 'S  3' precedes 'S  7'      
                            to_remove_10_2 = [is_S2(1:end-1) & is_S10(2:end), false];          % When 'S  2' precedes 'S  10'     
                            to_remove_8 = [is_S8(1:end-2), false, false];                      % When 'S  8' comes before the end  
                            to_remove_8R = is_8R;                                              % Remove TI triggers. 
                            to_remove = to_remove_1 | to_remove_2_3 | ...
                                        to_remove_3_4| to_remove_5_4 | ...
                                        to_remove_6_4 | to_remove_7_4 | ...
                                        to_remove_7_6 | to_remove_2_6 | ...
                                        to_remove_3_7 | to_remove_10_2 | ...
                                        to_remove_8 | to_remove_8R;
                            EEG.event(to_remove) = [];
                            disp([num2str(sum(to_remove)), ' spurious markers removed']);
                    
                            events = [EEG.event.latency];                                          
                            events = 1000*(events/EEG.srate);                                  % Convert event timings to ms                            
                            threshold = 0;                                                    
                            event_diff = [threshold, diff(events)];                            
                            reasonably_delayed_markers = event_diff >= threshold;              
                            markers_to_remove = event_diff < threshold;  
                            EEG.event(markers_to_remove) = [];
                            disp([num2str(sum(markers_to_remove)), ...                         
                               ' markers removed due to following the previous marker by <' ...
                               , num2str(threshold), 'ms'])
                            
                            if verb                                                            % Print remaining markers, latencies, and successive time differences
                                events = [EEG.event.latency];  
                                types = [EEG.event.type];
                                type_split = strsplit(types, ' ');     
                                events = events/1000;
                                combined = [num2cell(events)', type_split'];                   %#ok<NASGU> 
                                event_diff = [0, diff(events)];
                                combined = [num2cell(events)', type_split', ...
                                    num2cell(event_diff)'];
                                disp(combined);
                            end 
                        end 
                    
                        EEG = pop_resample(EEG, downsample_freq);                              % Downsample after marker extraction if non-cognitron task
                    end
                    
                    
                    %% Filter (DFT then high pass):
                    
                    ft_data = eeglab2fieldtrip(EEG, 'preprocessing');                          % DFT
                    cfg = [];
                    cfg.dftfilter = 'yes';
                    cfg.dftfreq = [50, 100, 150];
                    cfg.dftbandwidth = [1, 2, 3];
                    cfg.dftneighbourwidth = [2, 4, 6];
                    cfg.dftreplace = 'neighbour_fft';
                    ft_data = ft_preprocessing(cfg, ft_data);
                    
                    EEG.data = single(ft_data.trial{1});                                       % Replace EEG.data with the filtered data from FieldTrip
                    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'plotfreqz',0);                     % High pass filter: < 1Hz should be filtered out for ERP analyses
                    
                    
                    %% FFT (and IMD check): 
                    
                    if check_fft                   
                        n = length(EEG.data);
                        for i = 1:n_chan
                            yy = EEG.data(i, :);                                               % Extract time domain signal
                            f = (0:n/2)*(EEG.srate/n);                                         % Compute frequency vector
                            w = kaiser(n,2);                                                   % Create Kaiser filter
                            yy = yy(:).*w(:);                                                  % Apply Kaiser filter to signal
                            fft_values = fft(yy);                                              % Compute fft
                            fft_amp = abs(fft_values/n);                                       % Extract spectral amplitude (removing complex numbers)
                            fft_amp = fft_amp(1:n/2+1);                                        % Limit to positive frequencies only
                            fft_amp(2:end-1) = 2*fft_amp(2:end-1);                             % Compensate for removal of negative frequencies by doubling amplitude
                            figure('Position', [100, 100, 1000, 600]);                         % Plot
                            plot(f, fft_amp)
                            title([day, ', Session ', session, ', ', task, run_no, ', '  ...
                                'Channel ', EEG.chanlocs(i).labels, ' FFT']);
                            xlabel('Frequency (Hz)');
                            ylabel('Magnitude');
                            xlim([0 15]);                                                      % Interested particularly in 5Hz (IMD) range
                        end 
                    end
                    
                    %% Check number of trials:
                    
                    task_trial_no = containers.Map;                                            % Define expected trial numbers
                    task_trial_no('mst') = 60;
                    task_trial_no('fnat') = 18;
                    task_trial_no('f') = 27;
                    task_trial_no('wm') = 24;
    
                    task_trial_no('wwa') = 15;
                    task_trial_no('om') = 0;
                    task_trial_no('osa') = 0;
                    task_trial_no('pal') = 0;
                    task_trial_no('rs') = 0;
                    % Need to populate this for all tasks... 
                    
                    tmp_lab = {EEG.event.type};                                                % Trials may not be split evenly across EEG blocks   
                    n_trials = sum(strcmp(tmp_lab, 'S  2'));
                    if n_trials < task_trial_no(task)
                        error("error: fewer trials than expected")
                    elseif n_trials == task_trial_no(task)
                        disp("trial numbers as expected for task")
                    end 
                    
                    %% Save preprocessed file:

                    save(full_path, 'EEG');        
    
                end
    
            end
    
        end
    
    end

end