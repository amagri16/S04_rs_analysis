function [results,dscrp,idx] = find_PTE(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
start_index = size(results,2)+1;

sig1 = signal{1,1};
sig1_ch = labels{1,1};
sig2 = signal{2,1};
sig2_ch = labels{2,1};
sig3 = signal{3,1};
sig3_ch = labels{3,1};

max_delay = floor(0.8*split_len/10);

[~,eLag1,eDim1] = phaseSpaceReconstruction(sig1,'MaxDim',10,'MaxLag',max_delay);
[~,eLag2,eDim2] = phaseSpaceReconstruction(sig2,'MaxDim',10,'MaxLag',max_delay);
[~,eLag3,eDim3] = phaseSpaceReconstruction(sig3,'MaxDim',10,'MaxLag',max_delay);

DIM_S1 = eDim1;        % embedding dimension for the signal 1
DIM_S2 = eDim2;      % embedding dimension for the signal 2
DIM_S3 = eDim3;     % embedding dimension for the signal 3
STEP1 = eLag1;            % embedding delay
STEP2 = eLag2;            % embedding delay
STEP3 = eLag3;            % embedding delay
K = 4;              % number of nearest neighbors
WINDOW_RADIUS = 15;   % time-window radius for the partial TE estimates
LAGS = 0;         % lag values
FILTER_ORDER = 20;   % order of the post-processing moving average filter

% Function handle to the transfer entropy estimator
estimator = @(x, lag) transfer_entropy_pt(x(1,:), x(2,:), x(3,:), x(4,:), ...
    WINDOW_RADIUS, 'yLag', lag, 'k', K);
estimator_smooth = @(x, lag) filter((1/FILTER_ORDER)*ones(1,FILTER_ORDER),1, ...
    estimator(x, lag),[],2);


% feature description
dscrp{end+1} = string(sig1_ch)+'_'+string(sig3_ch)+'_to_'+string(sig2_ch)+'_Partial_Transfer_Entropy';
dscrp{end+1} = string(sig2_ch)+'_'+string(sig3_ch)+'_to_'+string(sig1_ch)+'_Partial_Transfer_Entropy';
dscrp{end+1} = string(sig1_ch)+'_'+string(sig2_ch)+'_to_'+string(sig3_ch)+'_Partial_Transfer_Entropy';
% feature calculation container for each paired-channels
TE1_all = [];
TE2_all = [];
TE3_all = [];
% epoching the signal
sig1_split = {};
sig2_split = {};
sig3_split = {};
num_sig1 = floor(size(sig1,1)/split_len);
num_sig2 = floor(size(sig2,1)/split_len);
num_sig3 = floor(size(sig3,1)/split_len);
% split the data into sub parts
for i=0:1:num_sig1-1
    sig1_split{1,i+1}=sig1(i*split_len+1:1:i*split_len+(split_len))';
end

for i=0:1:num_sig2-1
    sig2_split{1,i+1}=sig2(i*split_len+1:1:i*split_len+(split_len))';
end

for i=0:1:num_sig3-1
    sig3_split{1,i+1}=sig3(i*split_len+1:1:i*split_len+(split_len))';
end

% Assess information flow in the direction  1 <- 2
S1 = delay_embed(sig2_split, DIM_S1, STEP1);
S2 = delay_embed(sig1_split, DIM_S2, STEP2);
S3 = delay_embed(sig1_split, DIM_S3, STEP3);
W1 = delay_embed_future(S1, STEP1);
W2 = delay_embed_future(S2, STEP2);
W3 = delay_embed_future(S3, STEP3);

te1 = estimator_smooth([S1;S2;S3;W1],LAGS);
TE1_all = [TE1_all nanmean(sum(te1,1))];
te2 = estimator_smooth([S2;S1;S3;W2],LAGS);
TE2_all = [TE2_all nanmean(sum(te2,1))];
te3 = estimator_smooth([S3;S1;S2;W3],LAGS);
TE3_all = [TE3_all nanmean(sum(te3,1))];


results = [results mean(TE2_all)];
results = [results mean(TE1_all)];
results = [results mean(TE3_all)];


end_index = size(results,2);
idx = [start_index,end_index];