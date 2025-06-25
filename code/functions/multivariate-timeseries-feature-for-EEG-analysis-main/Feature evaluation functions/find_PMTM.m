function [results,dscrp,idx] = find_PMTM(y,f_s,ch_labels,split_length,results,dscrp)
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
beta_bd_1 = [12,18];
beta_bd_2 = [18,24];
beta_bd_3 = [24,30];

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman_pmtm';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_Spearman_pmtm';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_Spearman_pmtm';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_Spearman_pmtm';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_pmtm';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_pmtm_1';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_pmtm_2';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_pmtm_3';

        
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_delta_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_theta_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_alpha_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_beta_pmtm_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_theta_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_alpha_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_beta_pmtm_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_alpha_pmtm_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_beta_pmtm_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_beta_pmtm_power';

        % feature calculation container for each paired-channels
        % overall spearman's rank
        % Spearman_pmtm_all = [];
        % % spearman's rank on pmtm
        % delta_Spearman_pmtm_all = [];
        % theta_Spearman_pmtm_all = [];
        % alpha_Spearman_pmtm_all = [];
        % beta_Spearman_pmtm_all = [];
        % beta1_Spearman_pmtm_all = [];
        % beta2_Spearman_pmtm_all = [];
        % beta3_Spearman_pmtm_all = [];

        delta_delta_pmtm_power_all = [];
        delta_theta_pmtm_power_all = [];
        delta_alpha_pmtm_power_all = [];
        delta_beta_pmtm_power_all = [];

        theta_theta_pmtm_power_all = [];
        theta_alpha_pmtm_power_all = [];
        theta_beta_pmtm_power_all = [];

        alpha_alpha_pmtm_power_all = [];
        alpha_beta_pmtm_power_all = [];
        beta_beta_pmtm_power_all = []; 

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

            [pmtm1,freqs] = pmtm(sig1_diff,4,split_len,f_s);
            pmtm2 = pmtm(sig2_diff,4,split_len,f_s);
            
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            beta1_idx = freqs >= beta_bd_1(1) & freqs <= beta_bd_1(2);
            beta2_idx = freqs >= beta_bd_2(1) & freqs <= beta_bd_2(2);
            beta3_idx = freqs >= beta_bd_3(1) & freqs <= beta_bd_3(2);

            pmtm1 = pmtm1/sum(pmtm1);
            pmtm2 = pmtm2/sum(pmtm2);

            % Compute the spectral flatness
            % compute pearson
            % Spearman_pmtm = corr(pmtm1,pmtm2, 'Type', 'Spearman');
            % Spearman_pmtm_all = [Spearman_pmtm_all Spearman_pmtm];
            % 
            % delta_Spearman_pmtm = corr(pmtm1(delta_idx),pmtm2(delta_idx), 'Type', 'Spearman');
            % theta_Spearman_pmtm = corr(pmtm1(theta_idx),pmtm2(theta_idx), 'Type', 'Spearman');
            % alpha_Spearman_pmtm = corr(pmtm1(alpha_idx),pmtm2(alpha_idx), 'Type', 'Spearman');          
            % beta_Spearman_pmtm = corr(pmtm1(beta_idx),pmtm2(beta_idx), 'Type', 'Spearman');
            % beta1_Spearman_pmtm = corr(pmtm1(beta1_idx),pmtm2(beta1_idx), 'Type', 'Spearman');
            % beta2_Spearman_pmtm = corr(pmtm1(beta2_idx),pmtm2(beta2_idx), 'Type', 'Spearman');
            % beta3_Spearman_pmtm = corr(pmtm1(beta3_idx),pmtm2(beta3_idx), 'Type', 'Spearman');
            % 
            % delta_Spearman_pmtm_all = [delta_Spearman_pmtm_all abs(delta_Spearman_pmtm)];
            % theta_Spearman_pmtm_all = [theta_Spearman_pmtm_all abs(theta_Spearman_pmtm)];
            % alpha_Spearman_pmtm_all = [alpha_Spearman_pmtm_all abs(alpha_Spearman_pmtm)];
            % beta_Spearman_pmtm_all = [beta_Spearman_pmtm_all abs(beta_Spearman_pmtm)];
            % beta1_Spearman_pmtm_all = [beta1_Spearman_pmtm_all abs(beta1_Spearman_pmtm)];
            % beta2_Spearman_pmtm_all = [beta2_Spearman_pmtm_all abs(beta2_Spearman_pmtm)];
            % beta3_Spearman_pmtm_all = [beta3_Spearman_pmtm_all abs(beta3_Spearman_pmtm)]; 
            
            % power
            delta_power_sig1 = mean(pmtm1(delta_idx));
            theta_power_sig1 = mean(pmtm1(theta_idx));
            alpha_power_sig1 = mean(pmtm1(alpha_idx));
            beta_power_sig1 = mean(pmtm1(beta_idx));

            delta_power_sig2 = mean(pmtm2(delta_idx));
            theta_power_sig2 = mean(pmtm2(theta_idx));
            alpha_power_sig2 = mean(pmtm2(alpha_idx));
            beta_power_sig2 = mean(pmtm2(beta_idx));

            delta_delta_pmtm_power_all = [delta_delta_pmtm_power_all log(delta_power_sig1/delta_power_sig2)];
            delta_theta_pmtm_power_all = [delta_theta_pmtm_power_all log(delta_power_sig1/theta_power_sig2)];
            delta_alpha_pmtm_power_all = [delta_alpha_pmtm_power_all log(delta_power_sig1/alpha_power_sig2)];
            delta_beta_pmtm_power_all = [delta_beta_pmtm_power_all log(delta_power_sig1/beta_power_sig2)];

            theta_theta_pmtm_power_all = [theta_theta_pmtm_power_all log(theta_power_sig1/theta_power_sig2)];
            theta_alpha_pmtm_power_all = [theta_alpha_pmtm_power_all log(theta_power_sig1/alpha_power_sig2)];
            theta_beta_pmtm_power_all = [theta_beta_pmtm_power_all log(theta_power_sig1/beta_power_sig2)];

            alpha_alpha_pmtm_power_all = [alpha_alpha_pmtm_power_all log(alpha_power_sig1/alpha_power_sig2)];
            alpha_beta_pmtm_power_all = [alpha_beta_pmtm_power_all log(alpha_power_sig1/beta_power_sig2)];

            beta_beta_pmtm_power_all = [beta_beta_pmtm_power_all log(beta_power_sig1/beta_power_sig2)];
        end
        % add features into the results
        
        % results = [results mean(Spearman_pmtm_all)];
        % results = [results mean(delta_Spearman_pmtm_all)];
        % results = [results mean(theta_Spearman_pmtm_all)];
        % results = [results mean(alpha_Spearman_pmtm_all)];
        % results = [results mean(beta_Spearman_pmtm_all)];
        % results = [results mean(beta1_Spearman_pmtm_all)];
        % results = [results mean(beta2_Spearman_pmtm_all)];
        % results = [results mean(beta3_Spearman_pmtm_all)];


        results = [results mean((delta_delta_pmtm_power_all))];
        results = [results mean((delta_theta_pmtm_power_all))];
        results = [results mean((delta_alpha_pmtm_power_all))];
        results = [results mean((delta_beta_pmtm_power_all))];

        results = [results mean((theta_theta_pmtm_power_all))];
        results = [results mean((theta_alpha_pmtm_power_all))];
        results = [results mean((theta_beta_pmtm_power_all))];

        results = [results mean((alpha_alpha_pmtm_power_all))];
        results = [results mean((alpha_beta_pmtm_power_all))];

        results = [results mean((beta_beta_pmtm_power_all))];

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];