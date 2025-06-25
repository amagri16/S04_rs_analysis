function [results, dscrp, idx] = find_SL_BP(y,channel_labels, split_len,filters,results,dscrp,band_name)
 
signal = cell2mat(y');
CH = size(signal, 2); % Number of channels
EEG_chs = channel_labels;

signal_band = [];
for k = 1:CH-1
    signal_bd_now = filtfilt(filters, signal(:,k));
    signal_band = [signal_band signal_bd_now];
end
signal_band = [signal_band signal(:,CH)];

% Split the signals into segments
signal_segments = split_signal(signal_band, split_len);

% Initialize temporary results
temp_results = [];
start_index = size(results,2)+1;

for seg_idx = 1:length(signal_segments)
    % Get the segment signal
    seg_signal = signal_segments{seg_idx};

    % Calculate synchronization likelihood for the segment
    [S_matrix, S] = synchronization_likelyhood(seg_signal);
    
    % Store segment results in temporary results array
    seg_results = [];
    for i = 1:CH
        for j = i:CH
            if i ~= j
                seg_results = [seg_results, S_matrix(i, j)];
            end
        end
    end
    %seg_results = [seg_results, S];
    
    temp_results = [temp_results; seg_results];
end

% Average the computed features across all segments
avg_results = mean(temp_results, 1);

% Append the averaged features and their descriptions to the results and dscrp
results = [results, avg_results];
for i = 1:CH
    for j = i:CH
        if i ~= j
            dscrp{end+1} = [EEG_chs{i}, '_vs_', EEG_chs{j}, '_synchronization_' band_name];
        end
    end
end
%dscrp{end+1} = 'global_synchronization';

end_index = size(results,2);
idx = [start_index,end_index];
function segments = split_signal(signal, split_len)
    % Split the input signal into segments of length split_len
    n_segments = floor(size(signal, 1) / split_len);
    segments = cell(1, n_segments);

    for i = 1:n_segments
        start_idx = (i-1) * split_len + 1;
        end_idx = start_idx + split_len - 1;
        segments{i} = signal(start_idx:end_idx, :);
    end
end

end
