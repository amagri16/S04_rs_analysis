function [results,dscrp,idx] = find_TTE(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
results = results;
dscrp = dscrp;
start_index = size(results,2)+1;


%[~,eLag,eDim] = phaseSpaceReconstruction(cell2mat(y'));

DIM_S1 = 3;        % embedding dimension for the signal 1
DIM_S2 = 3;      % embedding dimension for the signal 2
DIM_FUTURE = 3;      % embed. dim. for the future of the signal
STEP = 1;            % embedding delay
K = 2;              % number of nearest neighbors
WINDOW_RADIUS = 5;   % time-window radius for the temporal TE estimates
LAGS = 1:5;         % lag values
FILTER_ORDER = 20;   % order of the post-processing moving average filter

% Function handle to the transfer entropy estimator
estimator = @(x) transfer_entropy_t(x(1,:), x(2,:), x(3,:), ...
    WINDOW_RADIUS, 'yLag', LAGS, 'k', K);
estimator_smooth = @(x) filter((1/FILTER_ORDER)*ones(1,FILTER_ORDER),1, ...
    estimator(x),[],2);



for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Temporal_Transfer_Entropy';
        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Temporal_Transfer_Entropy';
        % feature calculation container for each paired-channels
        TE12_all = [];
        TE21_all = [];
        % epoching the signal
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1,1)/split_len);
        num_sig2 = floor(size(sig2,1)/split_len);
        % split the data into sub parts
        for i=0:1:num_sig1-1
            sig1_split{1,i+1}=sig1(i*split_len+1:1:i*split_len+(split_len))';
        end

        for i=0:1:num_sig2-1
            sig2_split{1,i+1}=sig2(i*split_len+1:1:i*split_len+(split_len))';
        end

        % Assess information flow in the direction  2 <- 1
        S1 = delay_embed(sig2_split, DIM_S1, STEP);
        S2 = delay_embed(sig1_split, DIM_S2, STEP);
        W = delay_embed_future(S1, STEP);
        te21 = estimator_smooth([S1;S2;W]);
        TE21_all = [TE21_all nanmean(sum(te21,1))];

        % Assess information flow in the direction  1 <- 2
        S1 = delay_embed(sig1_split, DIM_S1, STEP);
        S2 = delay_embed(sig2_split, DIM_S2, STEP);
        W = delay_embed_future(S1, STEP);
        te12 = estimator_smooth([S1;S2;W]);
        TE12_all = [TE12_all nanmean(sum(te12,1))];

        results = [results mean(TE21_all)];
        results = [results mean(TE12_all)];

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];