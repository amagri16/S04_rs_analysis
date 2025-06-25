function [results,dscrp,idx] = find_PSDs(y,ch_labels,split_length,results,dscrp)
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
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_spectral_flatness';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_Spearman_psd_3';

        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_psd';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_psd';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_psd';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_1';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_2';
        % dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_psd_3';
        % feature calculation container for each paired-channels
        % overall spearman's rank
        Spearman_psd_all = [];
        % spectral flatness
        % spectral_flatness_sig1 = [];
        % spectral_flatness_sig2 = [];
        % 
        % delta_psd1_all = [];
        % theta_psd1_all = [];
        % alpha_psd1_all = [];
        % beta_psd1_all = [];
        % beta1_psd1_all = [];
        % beta2_psd1_all = [];
        % beta3_psd1_all = [];
        % 
        % delta_psd2_all = [];
        % theta_psd2_all = [];
        % alpha_psd2_all = [];
        % beta_psd2_all = [];
        % beta1_psd2_all = [];
        % beta2_psd2_all = [];
        % beta3_psd2_all = [];
        % spearman's rank on psd
        delta_Spearman_psd_all = [];
        theta_Spearman_psd_all = [];
        alpha_Spearman_psd_all = [];
        beta_Spearman_psd_all = [];
        beta1_Spearman_psd_all = [];
        beta2_Spearman_psd_all = [];
        beta3_Spearman_psd_all = [];
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
            [psd1,~] = pwelch(sig1_diff, [], [], 2560, 256);
            [psd2,~] = pwelch(sig2_diff, [], [], 2560, 256);
            
            psd1 = psd1/sum(psd1);
            psd2 = psd2/sum(psd2);


            % Compute the spectral flatness
            % spectral_flatness_val_sig1 = geomean(psd1) / mean(psd1);
            % spectral_flatness_val_sig2 = geomean(psd2) / mean(psd2);
            % 
            % spectral_flatness_sig1 = [spectral_flatness_sig1; spectral_flatness_val_sig1];
            % spectral_flatness_sig2 = [spectral_flatness_sig2; spectral_flatness_val_sig2];

            % compute pearson
            Spearman_psd = corr(psd1,psd2, 'Type', 'Spearman');
            Spearman_psd_all = [Spearman_psd_all Spearman_psd];
            
            delta_Spearman_psd = corr(psd1([11:1:41]),psd2([11:1:41]), 'Type', 'Spearman');
            theta_Spearman_psd = corr(psd1([41:1:81]),psd2([41:1:81]), 'Type', 'Spearman');
            alpha_Spearman_psd = corr(psd1([81:1:121]),psd2([81:1:121]), 'Type', 'Spearman');
            beta_Spearman_psd = corr(psd1([121:1:301]),psd2([121:1:301]), 'Type', 'Spearman');
            beta1_Spearman_psd = corr(psd1([121:1:181]),psd2([121:1:181]), 'Type', 'Spearman');
            beta2_Spearman_psd = corr(psd1([181:1:241]),psd2([181:1:241]), 'Type', 'Spearman');
            beta3_Spearman_psd = corr(psd1([241:1:301]),psd2([241:1:301]), 'Type', 'Spearman');

            delta_Spearman_psd_all = [delta_Spearman_psd_all abs(delta_Spearman_psd)];
            theta_Spearman_psd_all = [theta_Spearman_psd_all abs(theta_Spearman_psd)];
            alpha_Spearman_psd_all = [alpha_Spearman_psd_all abs(alpha_Spearman_psd)];
            beta_Spearman_psd_all = [beta_Spearman_psd_all abs(beta_Spearman_psd)];
            beta1_Spearman_psd_all = [beta1_Spearman_psd_all abs(beta1_Spearman_psd)];
            beta2_Spearman_psd_all = [beta2_Spearman_psd_all abs(beta2_Spearman_psd)];
            beta3_Spearman_psd_all = [beta3_Spearman_psd_all abs(beta3_Spearman_psd)]; 
            % % log on ratio
            % delta_psd1_all = [delta_psd1_all; log(mean(psd1([11:1:41])))];
            % theta_psd1_all = [theta_psd1_all; log(mean(psd1([41:1:81])))];
            % alpha_psd1_all = [alpha_psd1_all; log(mean(psd1([81:1:121])))];
            % beta_psd1_all = [beta_psd1_all; log(mean(psd1([121:1:301])))];
            % beta1_psd1_all = [beta1_psd1_all; log(mean(psd1([121:1:181])))];
            % beta2_psd1_all = [beta2_psd1_all;log(mean(psd1([181:1:241])))];
            % beta3_psd1_all = [beta3_psd1_all; log(mean(psd1([241:1:301])))];
            % 
            % delta_psd2_all = [delta_psd2_all; log(mean(psd2([11:1:41])))];
            % theta_psd2_all = [theta_psd2_all; log(mean(psd2([41:1:81])))];
            % alpha_psd2_all = [alpha_psd2_all; log(mean(psd2([81:1:121])))];
            % beta_psd2_all = [beta_psd2_all; log(mean(psd2([121:1:301])))];
            % beta1_psd2_all = [beta1_psd2_all; log(mean(psd2([121:1:181])))];
            % beta2_psd2_all = [beta2_psd2_all; log(mean(psd2([181:1:241])))];
            % beta3_psd2_all = [beta3_psd2_all; log(mean(psd2([241:1:301])))];
        end
        % add features into the results
        % p_spectral_flatness = corr(spectral_flatness_sig1, spectral_flatness_sig2,'Type', 'Spearman');
        % results = [results p_spectral_flatness];

        results = [results mean(Spearman_psd_all)];
        results = [results mean(delta_Spearman_psd_all)];
        results = [results mean(theta_Spearman_psd_all)];
        results = [results mean(alpha_Spearman_psd_all)];
        results = [results mean(beta_Spearman_psd_all)];
        results = [results mean(beta1_Spearman_psd_all)];
        results = [results mean(beta2_Spearman_psd_all)];
        results = [results mean(beta3_Spearman_psd_all)];

        % delta_psd = corr(delta_psd1_all,delta_psd2_all,'Type', 'Spearman');
        % theta_psd = corr(theta_psd1_all,theta_psd2_all,'Type', 'Spearman');
        % alpha_psd = corr(alpha_psd1_all,alpha_psd2_all,'Type', 'Spearman');
        % beta_psd = corr(beta_psd1_all,beta_psd2_all,'Type', 'Spearman');
        % beta1_psd = corr(beta1_psd1_all,beta1_psd2_all,'Type', 'Spearman');
        % beta2_psd = corr(beta2_psd1_all,beta2_psd2_all,'Type', 'Spearman');
        % beta3_psd = corr(beta3_psd1_all,beta3_psd2_all,'Type', 'Spearman');
        % 
        % results = [results delta_psd];
        % results = [results theta_psd];
        % results = [results alpha_psd];
        % results = [results beta_psd];
        % results = [results beta1_psd];
        % results = [results beta2_psd];
        % results = [results beta3_psd];

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];