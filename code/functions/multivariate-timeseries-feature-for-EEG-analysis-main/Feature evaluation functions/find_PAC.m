function [results, dscrp, idx] = find_PAC(y_cell,filters,channel_labels,split_len,results,dscrp)

% Input:
% y_cell: cell array containing signals for each channel
% Fs: the sampling frequency
% split_len: length of the signal segments for splitting
EEG_chs = channel_labels;
% Initialize results and dscrp
start_index = size(results,2)+1;
% Loop through all combinations of channels
n_channels = length(y_cell);
for ch1 = 1:n_channels
    for ch2 = 1:n_channels
        if ch1 == ch2
            continue;
        end

        % Split the signals into segments
        y1_segments = split_signal(y_cell{ch1}, split_len);
        y2_segments = split_signal(y_cell{ch2}, split_len);

        % Initialize temporary results
        temp_results = [];
        
        % Compute the PAC features for each segment and store in temp_results
        for seg_idx = 1:length(y1_segments)
            [ft_ch, ft_dscrp_ch] = compute_PAC(y1_segments{seg_idx}, y2_segments{seg_idx}, filters);
            temp_results = [temp_results; ft_ch];
        end
        
        % Average the computed features across all segments
        avg_ft_ch = mean(temp_results, 1);

        % Append the averaged features and their descriptions to the results and dscrp
        results = [results, avg_ft_ch];
        dscrp_ch_prefixed = strcat(ft_dscrp_ch,sprintf('_P:_%s_A:_%s', EEG_chs{ch1}, EEG_chs{ch2}));
        dscrp = [dscrp, dscrp_ch_prefixed];
    end
end
end_index = size(results,2);
idx = [start_index,end_index];

function segments = split_signal(signal, split_len)
    % Split the input signal into segments of length split_len
    n_segments = floor(length(signal) / split_len);
    segments = cell(1, n_segments);
    
    for i = 1:n_segments
        start_idx = (i-1) * split_len + 1;
        end_idx = start_idx + split_len - 1;
        segments{i} = signal(start_idx:end_idx);
    end
end
end