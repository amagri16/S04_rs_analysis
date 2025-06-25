function [results,dscrp,idx] = find_FFT(y,ch_labels,split_length,results,dscrp)
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
for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_alpha_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_fft';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_fft_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_fft_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_fft_3';

        
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_delta_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_theta_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_alpha_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_beta_fft_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_theta_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_alpha_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_beta_fft_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_alpha_fft_power';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_beta_fft_power';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_beta_fft_power';
        
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_alpha_combined_coupling_fft_power';


        % feature calculation container for each paired-channels
        % overall spearman's rank
        Spearman_fft_all = [];
        % spearman's rank on fft
        delta_Spearman_fft_all = [];
        theta_Spearman_fft_all = [];
        alpha_Spearman_fft_all = [];
        theta_alpha_Spearman_fft_all = [];
        beta_Spearman_fft_all = [];
        beta1_Spearman_fft_all = [];
        beta2_Spearman_fft_all = [];
        beta3_Spearman_fft_all = [];

        delta_delta_fft_power_all = [];
        delta_theta_fft_power_all = [];
        delta_alpha_fft_power_all = [];
        delta_beta_fft_power_all = [];

        theta_theta_fft_power_all = [];
        theta_alpha_fft_power_all = [];
        theta_beta_fft_power_all = [];

        alpha_alpha_fft_power_all = [];
        alpha_beta_fft_power_all = [];
        beta_beta_fft_power_all = []; 

        theta_alpha_combined_fft_power_all = [];

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
            fft1 = fft(sig1_diff,2560);
            fft2 = fft(sig2_diff,2560);

            fft1 = abs(fft1(1:1280+1));
            fft2 = abs(fft2(1:1280+1));

            fft1(2:end-1) = 2*abs(fft1(2:end-1));
            fft2(2:end-1) = 2*abs(fft2(2:end-1));

            fft1 = fft1/sum(fft1);
            fft2 = fft2/sum(fft2);

            % Compute the spectral flatness
            % compute pearson
            Spearman_fft = corr(fft1,fft2, 'Type', 'Spearman');
            Spearman_fft_all = [Spearman_fft_all Spearman_fft];
            
            delta_Spearman_fft = corr(fft1([11:1:41]),fft2([11:1:41]), 'Type', 'Spearman');
            theta_Spearman_fft = corr(fft1([41:1:81]),fft2([41:1:81]), 'Type', 'Spearman');
            alpha_Spearman_fft = corr(fft1([81:1:121]),fft2([81:1:121]), 'Type', 'Spearman');          
            theta_alpha_Spearman_fft = corr(fft1([41:1:121]),fft2([41:1:121]), 'Type', 'Spearman');
            beta_Spearman_fft = corr(fft1([121:1:301]),fft2([121:1:301]), 'Type', 'Spearman');
            beta1_Spearman_fft = corr(fft1([121:1:181]),fft2([121:1:181]), 'Type', 'Spearman');
            beta2_Spearman_fft = corr(fft1([181:1:241]),fft2([181:1:241]), 'Type', 'Spearman');
            beta3_Spearman_fft = corr(fft1([241:1:301]),fft2([241:1:301]), 'Type', 'Spearman');

            delta_Spearman_fft_all = [delta_Spearman_fft_all abs(delta_Spearman_fft)];
            theta_Spearman_fft_all = [theta_Spearman_fft_all abs(theta_Spearman_fft)];
            alpha_Spearman_fft_all = [alpha_Spearman_fft_all abs(alpha_Spearman_fft)];
            theta_alpha_Spearman_fft_all = [theta_alpha_Spearman_fft_all abs(theta_alpha_Spearman_fft)];
            beta_Spearman_fft_all = [beta_Spearman_fft_all abs(beta_Spearman_fft)];
            beta1_Spearman_fft_all = [beta1_Spearman_fft_all abs(beta1_Spearman_fft)];
            beta2_Spearman_fft_all = [beta2_Spearman_fft_all abs(beta2_Spearman_fft)];
            beta3_Spearman_fft_all = [beta3_Spearman_fft_all abs(beta3_Spearman_fft)]; 
            
            % power
            delta_power_sig1 = mean(fft1([11:1:41]));
            theta_power_sig1 = mean(fft1([41:1:81]));
            alpha_power_sig1 = mean(fft1([81:1:121]));
            beta_power_sig1 = mean(fft1([121:1:301]));

            delta_power_sig2 = mean(fft2([11:1:41]));
            theta_power_sig2 = mean(fft2([41:1:81]));
            alpha_power_sig2 = mean(fft2([81:1:121]));
            beta_power_sig2 = mean(fft2([121:1:301]));

            theta_alpha_combined_sig1 = mean(fft1(41:1:121)); 
            theta_alpha_combined_sig2 = mean(fft2(41:1:121)); 

            delta_delta_fft_power_all = [delta_delta_fft_power_all log(delta_power_sig1/delta_power_sig2)];
            delta_theta_fft_power_all = [delta_theta_fft_power_all log(delta_power_sig1/theta_power_sig2)];
            delta_alpha_fft_power_all = [delta_alpha_fft_power_all log(delta_power_sig1/alpha_power_sig2)];
            delta_beta_fft_power_all = [delta_beta_fft_power_all log(delta_power_sig1/beta_power_sig2)];

            theta_theta_fft_power_all = [theta_theta_fft_power_all log(theta_power_sig1/theta_power_sig2)];
            theta_alpha_fft_power_all = [theta_alpha_fft_power_all log(theta_power_sig1/alpha_power_sig2)];
            theta_beta_fft_power_all = [theta_beta_fft_power_all log(theta_power_sig1/beta_power_sig2)];

            alpha_alpha_fft_power_all = [alpha_alpha_fft_power_all log(alpha_power_sig1/alpha_power_sig2)];
            alpha_beta_fft_power_all = [alpha_beta_fft_power_all log(alpha_power_sig1/beta_power_sig2)];

            beta_beta_fft_power_all = [beta_beta_fft_power_all log(beta_power_sig1/beta_power_sig2)];
            theta_alpha_combined_fft_power_all = [theta_alpha_combined_fft_power_all log(theta_alpha_combined_sig1/theta_alpha_combined_sig2)];
        end
        % add features into the results

        results = [results mean(Spearman_fft_all)];
        results = [results mean(delta_Spearman_fft_all)];
        results = [results mean(theta_Spearman_fft_all)];
        results = [results mean(alpha_Spearman_fft_all)];
        results = [results mean(theta_alpha_Spearman_fft_all)];
        results = [results mean(beta_Spearman_fft_all)];
        results = [results mean(beta1_Spearman_fft_all)];
        results = [results mean(beta2_Spearman_fft_all)];
        results = [results mean(beta3_Spearman_fft_all)];


        results = [results mean(abs(delta_delta_fft_power_all))];
        results = [results mean(abs(delta_theta_fft_power_all))];
        results = [results mean(abs(delta_alpha_fft_power_all))];
        results = [results mean(abs(delta_beta_fft_power_all))];

        results = [results mean(abs(theta_theta_fft_power_all))];
        results = [results mean(abs(theta_alpha_fft_power_all))];
        results = [results mean(abs(theta_beta_fft_power_all))];

        results = [results mean(abs(alpha_alpha_fft_power_all))];
        results = [results mean(abs(alpha_beta_fft_power_all))];

        results = [results mean(abs(beta_beta_fft_power_all))];
        results = [results mean(abs(theta_alpha_combined_fft_power_all))];


        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];