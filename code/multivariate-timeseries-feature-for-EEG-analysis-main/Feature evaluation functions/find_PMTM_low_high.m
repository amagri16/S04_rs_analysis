function [results,dscrp,idx] = find_PMTM_low_high(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
start_index = size(results,2)+1;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [12,30];

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_low_high';

        % feature calculation container for each paired-channels
        low_high_pmtm_power_all = [];

        % epoching the signal
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1,1)/split_len);
        num_sig2 = floor(size(sig2,1)/split_len);
        % split the data into sub parts
        for i=0:1:num_sig1-1
            sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
        end

        for i=0:1:num_sig2-1
            sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
        end

        sig1_split_size = size(sig1_split,1);
        sig2_split_size = size(sig2_split,1);
        min_sig_size = min(sig1_split_size,sig2_split_size);
        % traverse through the splited signals

        for k = 1:1:min_sig_size
            sig1_now = sig1_split{k};
            sig2_now = sig2_split{k};
            
            % Compute the power spectral density of the signal using Welch's method

            sig1_diff = sig1_now(2:end) - sig1_now(1:end-1);
            sig2_diff = sig2_now(2:end) - sig2_now(1:end-1);

            [pmtm1,freqs] = pmtm(sig1_diff,4,split_len,256);
            pmtm2 = pmtm(sig2_diff,4,split_len,256);
            
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);

            pmtm1 = pmtm1/sum(pmtm1);
            pmtm2 = pmtm2/sum(pmtm2);

            % power
            delta_power_sig1 = mean(pmtm1(delta_idx));
            theta_power_sig1 = mean(pmtm1(theta_idx));
            alpha_power_sig2 = mean(pmtm2(alpha_idx));
            beta_power_sig2 = mean(pmtm2(beta_idx));

            low_high_pmtm_power_all = [low_high_pmtm_power_all log((delta_power_sig1+theta_power_sig1)...
                                                            /(alpha_power_sig2+beta_power_sig2))];
        end
        % add features into the results
        results = [results mean((low_high_pmtm_power_all))];
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];