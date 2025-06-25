
function [encoding_periods] = epoch_FNAT(subject, week, day, session, run_no, epoch_size)

    % epoch_FNAT: A function to epoch preprocessed face-name association
    % task EEG data. 

    epoch_encoding_periods = true;

    %% Establish direc
    direc = '/Users/simonwilliamson/Desktop/PhD/EEG analysis/';

    %% Load in preprocessed data:
    
    path1 = [direc, 'preprocessed', filesep, 'subject_0', num2str(subject),...
        filesep, 'FNAT', filesep, 'week_', week, filesep, day, ...
        '_Session_', session, '_fnat', run_no, '.mat'];
    
    disp(path1);
    
    data = load(path1);
    EEG = data.EEG;
    n_chan = size(EEG.data, 1) - 1;

    %% Load in behavioural data:

    switch day                                                                 % Find session number
    case 'pre_day'
        x = 1;
    case 'stim_day3'
        x = 2;
    case 'post_day'
        x = 3;
    case 'follow_up1'
        x = 4;
    case 'follow_up3'
        x = 5;
    otherwise
            error('Invalid day.');
    end
    
    if week == '1'
        x = num2str(x);
    else
        x = num2str(x + 5);
    end
    
    
    direc2 = ['/Users/simonwilliamson/Desktop/PhD/BEH analysis/', ...          % Create filepath
                'beh_fnat_preprocessed/subject_', num2str(subject), ...
                '/week_', week, '/'];
    beh_file = dir(fullfile(direc2, [num2str(subject), '_S_', x, '_*.mat']));
    path2 = [direc2, beh_file.name];
    disp(path2);
    
    beh_data = load(path2);
    
    accuracy = [beh_data.StructTask.Accuracy];
    
    if length(accuracy) ~= 54        
        error('Incorrect number of trials in behavioural data')
    end 
    
    switch run_no                                                              % Get trials corresponding to current run
        case '1'
            accuracy = accuracy(1:18);
            beh_events = [beh_data.StructTask(1:18)];
        case '2'
            accuracy = accuracy(19:36);
            beh_events = [beh_data.StructTask(19:36)];
        case '3'
            accuracy = accuracy(37:54);
            beh_events = [beh_data.StructTask(37:54)];
    end 

    %% Sync behavioral markers with EEG signal:

    all_types = {EEG.event.type};                     
    index = find(strcmp(all_types, 'S  6'), 1);   
    
    sorted_beh_markers = sort([beh_events.RespTimeStart]);
    eeg_delay = EEG.event(index).latency/1000 - sorted_beh_markers(1);    
    beh_events_common = [beh_events.RespTimeStart] + eeg_delay; 

    % Do they look synced?
    z = ones(size(beh_events_common));
    x = ones(size(EEG.event));
    plot([EEG.event.latency], x, 'rx')
    hold on;
    plot(beh_events_common*1000, z, 'bo')
    
    % Sync all relevant markers:
    encoding_start = ([beh_events.TrialStart] + eeg_delay)*1000;
    RespTimeStart = ([beh_events.RespTimeStart] + eeg_delay)*1000;
    FacePresented = ([beh_events.FacePresented] + eeg_delay)*1000;

    
    %% Create epochs for encoding periods based on behavioural markers:
    
    % Determine size of epochs:
    if length(epoch_size) ~= 2
        error("error: incorrect number of values provided for epoch_size argument - please provide two values")
    end 
    epoch_start = epoch_size(1);
    epoch_end = epoch_size(2);
    
    if epoch_encoding_periods == true

        % Find encoding periods:
        encoding_end = encoding_start + 4*EEG.srate;
        epoch_length = (epoch_end - epoch_start)*EEG.srate;                                              
        n_epochs = 18;
        
        % Round latencies:
        encoding_start = round(encoding_start);
        encoding_end = round(encoding_end);

        % Take epochs around the encoding event as specified in epoch_size:
        epoch_to_extract_start = encoding_start - epoch_start*EEG.srate;
        epoch_to_extract_end = epoch_to_extract_start + epoch_length;
        
        % Create epochs:
        epochs = zeros(n_chan, epoch_length, n_epochs);
        for j = 1:n_epochs
            epochs(:, :, j) = EEG.data(1:n_chan, epoch_to_extract_start(j):epoch_to_extract_end(j)-1);
        end
                                            
        
        disp(['Encoding period epochs structure created: ',...
            num2str(size(epochs, 3)), ' epochs of duration ',...
            num2str(size(epochs, 2)/EEG.srate), ...
            ' seconds in ', num2str(size(epochs, 1)), ' channels']);

        encoding_periods = epochs;

    end 

end