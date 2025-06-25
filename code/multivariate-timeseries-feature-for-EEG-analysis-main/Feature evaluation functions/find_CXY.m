function [results,dscrp,idx] = find_CXY(y,f_s,ch_labels,split_length,results,dscrp)

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
results = results;
dscrp = dscrp;

delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [12,30];
beta_bd_1 = [12,18];
beta_bd_2 = [18,24];
beta_bd_3 = [24,30];
gamma_bd = [30,35];
start_index = size(results,2)+1;

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_delta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_theta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_alpha_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_1';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_2';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_beta_coherence_3';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_gamma_coherence';
        %dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_coherence_flatness_all_bands';
        % feature calculation container for each paired-channels
        mean_delta_csd_all = [];
        mean_theta_csd_all = [];
        mean_alpha_csd_all = [];
        mean_beta_csd_all = [];
        mean_beta1_csd_all = [];
        mean_beta2_csd_all = [];
        mean_beta3_csd_all = [];
        mean_gamma_csd_all = [];
        %coherence_flatness = [];
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
            %feature calculation
            
            %compute coherence density
            % sig1_now = sig1_now(2:end) - sig1_now(1:end-1);
            % sig2_now = sig2_now(2:end) - sig2_now(1:end-1);

            [csd, freqs] = mscohere(sig1_now, sig2_now,[],[], 0:0.05:35, f_s);
            
            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            beta1_idx = freqs >= beta_bd_1(1) & freqs <= beta_bd_1(2);
            beta2_idx = freqs >= beta_bd_2(1) & freqs <= beta_bd_2(2);
            beta3_idx = freqs >= beta_bd_3(1) & freqs <= beta_bd_3(2);
            gamma_idx = freqs >= gamma_bd(1) & freqs <= gamma_bd(2);

            csd_delta = sum(csd(delta_idx));
            csd_theta = sum(csd(theta_idx));
            csd_alpha = sum(csd(alpha_idx));
            csd_beta = sum(csd(beta_idx));
            csd_beta1 = sum(csd(beta1_idx));
            csd_beta2 = sum(csd(beta2_idx));
            csd_beta3 = sum(csd(beta3_idx));
            csd_gamma = sum(csd(gamma_idx));

            mean_delta_csd_all = [mean_delta_csd_all csd_delta];
            mean_theta_csd_all = [mean_theta_csd_all csd_theta];
            mean_alpha_csd_all = [mean_alpha_csd_all csd_alpha];
            mean_beta_csd_all = [mean_beta_csd_all csd_beta];
            mean_beta1_csd_all = [mean_beta1_csd_all csd_beta1];
            mean_beta2_csd_all = [mean_beta2_csd_all csd_beta2];
            mean_beta3_csd_all = [mean_beta3_csd_all csd_beta3];
            mean_gamma_csd_all = [mean_gamma_csd_all csd_gamma];
            %coherence_flatness = [coherence_flatness geomean(csd)/mean(csd)];
        end
        % add features into the results
        results = [results mean(mean_delta_csd_all)];
        results = [results mean(mean_theta_csd_all)];
        results = [results mean(mean_alpha_csd_all)];
        results = [results mean(mean_beta_csd_all)];
        results = [results mean(mean_beta1_csd_all)];
        results = [results mean(mean_beta2_csd_all)];
        results = [results mean(mean_beta3_csd_all)];
        results = [results mean(mean_gamma_csd_all)];
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];