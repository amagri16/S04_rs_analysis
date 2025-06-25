%% Automated EEG Preprocessing Script
% This script automates the preprocessing of all EEG recordings
% for subject S04 across weeks 1 and 2.
% It will process all rs raw data files.
% A comprehensive markdown file is available in the /code/docs folder, with a step-by-step guide to the preprocessing steps.

clear
clc

%% Initialize EEGLAB and FieldTrip
addpath '/Users/andreamagri/Desktop/thesis/S04_analysis/code/toolboxes/eeglab2025.0.0';
eeglab;
addpath '/Users/andreamagri/Desktop/thesis/S04_analysis/code/toolboxes/fieldtrip-20250523';
ft_defaults

%% Define Parameters
verb = true;                                                               % Print information to terminal
downsample_freq = 1000;                                                    % Decide frequency to downsample to (Hz)
n_chan = 14;                                                               % Number of EEG channels (excluding LED)
check_fft = false;                                                         % Plot FFTs of Kaiser-windowed signal for each channel
subject = '04';                                                            % Subject number

% Define all sessions to process
sessions = {
    % Week 1
    struct('week', '1', 'day', 'pre_day', 'session', '', 'has_session_folder', false),
    struct('week', '1', 'day', 'post_day', 'session', '', 'has_session_folder', false),
    struct('week', '1', 'day', 'stim_day3', 'session', '1', 'has_session_folder', true),
    struct('week', '1', 'day', 'follow_up1', 'session', '', 'has_session_folder', false),
    struct('week', '1', 'day', 'follow_up3', 'session', '', 'has_session_folder', false),
    % Week 2
    struct('week', '2', 'day', 'pre_day', 'session', '', 'has_session_folder', false),
    struct('week', '2', 'day', 'post_day', 'session', '', 'has_session_folder', false),
    struct('week', '2', 'day', 'stim_day3', 'session', '1', 'has_session_folder', true),
    struct('week', '2', 'day', 'follow_up1', 'session', '', 'has_session_folder', false),
    struct('week', '2', 'day', 'follow_up3', 'session', '', 'has_session_folder', false)
};

% Define recordings for each session
recordings = {'1', '2'};

% Create a log file to track processing
log_file = fullfile('/Users/andreamagri/Desktop/thesis/S04_analysis/data/processed', 'preprocessing_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, 'EEG Preprocessing Log\n');
fprintf(fid, 'Started: %s\n\n', char(datetime('now')));

%% Process each session and recording
for s = 1:length(sessions)
    current_session = sessions{s};
    
    for r = 1:length(recordings)
        run_no = recordings{r};
        
        % Construct file paths
        if current_session.has_session_folder
            rawFolder = fullfile('/Users/andreamagri/Desktop/thesis/S04_analysis/data/raw', ...
                ['week' current_session.week], current_session.day, ['session_' current_session.session]);
        else
            rawFolder = fullfile('/Users/andreamagri/Desktop/thesis/S04_analysis/data/raw', ...
                ['week' current_session.week], current_session.day);
        end
        
        % Check for all required files
        hdrFile = ['rs' run_no '.vhdr'];
        eegFile = ['rs' run_no '.eeg'];
        vmrkFile = ['rs' run_no '.vmrk'];
        
        hdrPath = fullfile(rawFolder, hdrFile);
        eegPath = fullfile(rawFolder, eegFile);
        vmrkPath = fullfile(rawFolder, vmrkFile);
        
        % Check if all required files exist
        if ~exist(hdrPath, 'file')
            warning('Header file not found: %s', hdrPath);
            fprintf(fid, 'SKIPPED: Header file not found: %s\n', hdrPath);
            continue;
        end
        if ~exist(eegPath, 'file')
            warning('EEG data file not found: %s', eegPath);
            fprintf(fid, 'SKIPPED: EEG data file not found: %s\n', eegPath);
            continue;
        end
        if ~exist(vmrkPath, 'file')
            warning('Marker file not found: %s', vmrkPath);
            fprintf(fid, 'SKIPPED: Marker file not found: %s\n', vmrkPath);
            continue;
        end
        
        % Create session string for display
        session_str = '';
        if current_session.has_session_folder
            session_str = ['/session_' current_session.session];
        end
        
        % Check if output file already exists
        output_dir = fullfile('/Users/andreamagri/Desktop/thesis/S04_analysis/data/processed', ...
            ['week_' current_session.week], current_session.day);
        f_name = sprintf('S%s_%s_wk%s_rs%s.mat', ...
            subject, current_session.day, current_session.week, run_no);
        full_save_path = fullfile(output_dir, f_name);
        
        if exist(full_save_path, 'file')
            warning('Output file already exists: %s', full_save_path);
            fprintf(fid, 'SKIPPED: Output file already exists: %s\n', full_save_path);
            continue;
        end
        
        fprintf('Processing: Week %s, %s%s, Recording %s\n', ...
            current_session.week, current_session.day, session_str, run_no);
        fprintf(fid, 'Processing: Week %s, %s%s, Recording %s\n', ...
            current_session.week, current_session.day, session_str, run_no);
        
        try
            % Start timing
            tic;
            start_time = datetime('now');
            fprintf('Started at: %s\n', char(start_time));
            fprintf(fid, 'Started at: %s\n', char(start_time));
            
            % Load data
            [EEG, header] = pop_loadbv(rawFolder, hdrFile);
            
            if size(EEG.data,1) ~= (n_chan + 1)
                warning('Number of channels loaded (%d) does not match expected (%d).', ...
                    size(EEG.data,1), n_chan+1);
                fprintf(fid, 'WARNING: Channel count mismatch in %s\n', hdrPath);
            end
            
            % Define task type (always resting state for this script)
            task = 'rs';
            disp("resting state recording confirmed")
            cognitron = false;
            
            % Low-pass filter
            EEG.data = EEG.data/10;                                                    % Undo amplification
            EEG = pop_eegfiltnew(EEG, 'hicutoff', downsample_freq/4, 'plotfreqz', 0);  % Set filter cutoff to half the nyquist freq
            
            % Downsample (for non-cognitron tasks, do this after marker extraction)
            EEG = pop_resample(EEG, downsample_freq);
            
            % DFT and High-pass Filter
            ft_data = eeglab2fieldtrip(EEG, 'preprocessing');
            cfg = [];
            cfg.dftfilter = 'yes';
            cfg.dftfreq = [50, 100, 150];
            cfg.dftbandwidth = [1, 2, 3];
            cfg.dftneighbourwidth = [2, 4, 6];
            cfg.dftreplace = 'neighbour_fft';
            ft_data = ft_preprocessing(cfg, ft_data);
            
            EEG.data = single(ft_data.trial{1});
            EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'plotfreqz', 0);
            
            % Optional FFT check
            if check_fft
                n = length(EEG.data);
                for i = 1:n_chan
                    yy = EEG.data(i, :);
                    f = (0:n/2)*(EEG.srate/n);
                    w = kaiser(n,2);
                    yy = yy(:).*w(:);
                    fft_values = fft(yy);
                    fft_amp = abs(fft_values/n);
                    fft_amp = fft_amp(1:n/2+1);
                    fft_amp(2:end-1) = 2*fft_amp(2:end-1);
                    figure('Position', [100, 100, 1000, 600]);
                    plot(f, fft_amp)
                    title([current_session.day, ', Session ', current_session.session, ', rs', run_no, ', ' ...
                        'Channel ', EEG.chanlocs(i).labels, ' FFT']);
                    xlabel('Frequency (Hz)');
                    ylabel('Magnitude');
                    xlim([0 15]);
                end
            end
            
            % Check if output directory exists
            if ~exist(output_dir, 'dir')
                error('Output directory does not exist: %s\nPlease create the directory structure first.', output_dir);
            end
            
            % Save preprocessed data
            save(full_save_path, 'EEG');
            
            % End timing
            end_time = datetime('now');
            processing_time = toc;
            fprintf('Completed at: %s (took %.2f seconds)\n', char(end_time), processing_time);
            fprintf(fid, 'Completed at: %s (took %.2f seconds)\n', char(end_time), processing_time);
            
            fprintf('Successfully processed and saved: %s\n', f_name);
            fprintf(fid, 'SUCCESS: Processed and saved: %s\n', f_name);
            
        catch ME
            % End timing even if there's an error
            end_time = datetime('now');
            processing_time = toc;
            fprintf('Failed at: %s (after %.2f seconds)\n', char(end_time), processing_time);
            fprintf(fid, 'Failed at: %s (after %.2f seconds)\n', char(end_time), processing_time);
            
            warning('Error processing %s: %s', fullfile(rawFolder, hdrFile), ME.message);
            fprintf(fid, 'ERROR: Failed to process %s: %s\n', hdrPath, ME.message);
            continue;
        end
    end
end

fprintf(fid, '\nPreprocessing completed: %s\n', char(datetime('now')));
fclose(fid);
fprintf('\nPreprocessing complete! Check preprocessing_log.txt for details.\n'); 