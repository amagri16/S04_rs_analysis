%% TIAD2 EEG Preprocessing


% Current sections
  % - Load in Data
  % - Define the task
  % - Low-pass filter
  % - Extract markers and downsample
  % - DFT and High-pass Filter
  % - FFT
  % - Check number of trials


%% Clear workspace
clear
clc


%% Initialise eeglab and fieldtrip

% Add paths to EEGLAB and FieldTrip toolboxes within the code/toolboxes folder
addpath '/Users/andreamagri/Desktop/thesis/S04_analysis/code/toolboxes/eeglab2025.0.0';
eeglab;
addpath '/Users/andreamagri/Desktop/thesis/S04_analysis/code/toolboxes/fieldtrip-20250523';
ft_defaults


%% Set Parameters


verb = true;                                                               % Print information to terminal
downsample_freq = 1000;                                                    % Decide frequency to downsample to (Hz)
n_chan = 14;                                                               % Number of EEG channels (excluding LED)
check_fft = false;                                                          % Plot FFTs of Kaiser-windowed signal for each channel


% *** IMPORTANT: Set these parameters for the file you want to preprocess ***
% These parameters are used for saving and task-specific logic
subject = '04';   % e.g., '04', '05'
week = '1';       % e.g., '1', '2', '3', '4'
day = 'post_day'; % e.g., 'pre_day', 'stim_day1', 'follow_up3'
session = '';     % e.g., '1', '2'. Leave empty if no session folder (like pre_day)
task = 'rs';      % e.g., 'rs', 'ssvep', 'wm', 'om', 'pal', 'wwa', 'osa', 'fnat', 'f', 'mst'
run_no = '2';     % e.g., '1', '2', or 'd' for delayed, or '' for no number (i.e. wm, om)




%% Load in data:

% Load data using hardcoded path from user's example
% This assumes the raw data file is directly accessible at this path.
%  —————————————— FIXED “Load in data” ——————————————
rawFolder = '/Users/andreamagri/Desktop/thesis/S04_analysis/data/raw/week1/post_day';
hdrFile   = 'rs2.vhdr';
[EEG, header] = pop_loadbv(rawFolder, hdrFile);

if size(EEG.data,1) ~= (n_chan + 1)
    warning('Number of channels loaded (%d) does not match expected (%d).', ...
            size(EEG.data,1), n_chan+1);
end


% Define the task type based on parameters
if matches(task, 'mst') || matches(task, 'fnat') || matches(task, 'f')   
    disp("non-cognitron task confirmed")
    cognitron = false;
elseif matches(task, 'wm') || matches(task, 'om') || matches(task, 'pal') || matches(task, 'wwa') || matches(task, 'osa')
    disp("cognitron task confirmed")
    cognitron = true;
elseif matches(task, 'rs')
    disp("resting state recording confirmed")
    cognitron = false;
elseif matches(task, 'ssvep')
    disp("ssvep recording confirmed")
    cognitron = false;
else
    error("error: task type not found")
end



%% Low-pass filter the data:

EEG.data = EEG.data/10;                                                    % Undo amplification 
EEG = pop_eegfiltnew(EEG, 'hicutoff',downsample_freq/4,'plotfreqz',0);     % Set filter cutoff to half the nyquist freq. of downsampled sampling rate


%% Extract markers and Downsample:

% Marker extraction logic based on provided script content
if cognitron == true


    EEG = pop_resample(EEG, downsample_freq);                              % Downsample first if cognitron task
    LED_signal = EEG.data(n_chan+1, :);                                    % Isolate the LED signal
    LED_deriv1 = [0, diff(LED_signal)];                                    % Find 1st derivative peaks and troughs
    threshold = mean(LED_deriv1) + 10*std(LED_deriv1);                     % Set a reasonable threshold for finding peaks
    if matches(task, 'wwa')
        threshold = mean(LED_deriv1) + 4*std(LED_deriv1);   
    end
    [peaks, peak_times] = findpeaks(LED_deriv1, ...
        'MinPeakHeight', threshold);
    LED_deriv1_inverted = -LED_deriv1;                                     % plot(LED_deriv1)  % plot(LED_deriv1_inverted)
    [troughs, trough_times] = findpeaks(LED_deriv1_inverted, ...
        'MinPeakHeight', threshold);

    if matches(task, 'wwa') & matches(run_no, '1')
        if length(peaks) > 16                                              % For WWA1, there should be 16 troughs (stimuli)
            peaks = peaks(1:16);                                           % and 16 peaks (responses)
            peak_times = peak_times(1:16);
        end
        if length(troughs) > 16
            troughs = troughs(1:16);
            trough_times = trough_times(1:16);
        end
    end
   
    if length(troughs) == length(peaks)-1                                  % Remove last peak if spurious
        peaks = peaks(1:end-1);
    end 
    if length(troughs) ~= length(peaks)
        error("error: The numbers of stimuli" + ...
            " and response markers are not equal")
    end 
    
    disp("Adding optical events to EEG data structure")                    % Add markers to EEG structure 
    for i = 1:length(trough_times)                                         % Add each trial in turn (peaks are responses, troughs are stimuli / trial starts)
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
        fft_amp = fft_amp(1:n/2+1);
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


%% Check number of trials (only need to run for tasks):


%task_trial_no = containers.Map;                                            % Define expected trial numbers
%task_trial_no('mst') = 60;
%task_trial_no('fnat') = 18;
%task_trial_no('f') = 27;
%task_trial_no('wm') = 24;
%task_trial_no('wwa') = 16;
% Need to populate this for all tasks... 


%tmp_lab = {EEG.event.type};                                                % Trials may not be split evenly across EEG blocks   
%n_trials = sum(strcmp(tmp_lab, 'S  2'));
%if n_trials < task_trial_no(task)
%    error("error: fewer trials than expected")
%else
%    disp("trial numbers as expected for task")
%end 


%% Save preprocessed file:

% Define the base directory for saving preprocessed data
project_root = '/Users/andreamagri/Desktop/thesis/S04_analysis';
preprocessed_data_base_direc = fullfile(project_root, 'data', 'processed');

% Define the folder path for saving preprocessed data based on parameters
% Assuming saving to data/processed/week_X/day_X/
output_dir = fullfile(preprocessed_data_base_direc, ['week_' week], day);

% Create the folder if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create the filename based on parameters
% Filename format: S_subject_day_wk_week_task_run_no.mat (e.g., S04_post_day_wk1_rs1.mat)
f_name = sprintf('S%s_%s_wk%s_%s%s.mat', subject, day, week, task, run_no);

% Combine the folder path with the file name
full_save_path = fullfile(output_dir, f_name);

% Save the EEG structure
save(full_save_path, 'EEG');

