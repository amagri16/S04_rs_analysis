function [results,dscrp,idx] = find_Cross_Bispectrum(y,Fs,ch_labels,split_length,results,dscrp)

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
results = results;
dscrp = dscrp;
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
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_delta';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_theta';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_alpha';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_beta'; 
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_delta_theta';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_delta_alpha';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_delta_beta';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_theta_alpha';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_theta_beta';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_alpha_beta';
    
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_delta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_theta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_alpha';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_beta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_delta_theta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_delta_alpha';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_delta_beta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_theta_alpha';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_theta_beta';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_CBS_alpha_beta';


%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_delta';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_theta';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_alpha';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_beta';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_beta1';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_beta2';
%         dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_CBS_beta3';

        % feature calculation container for each paired-channels
        CBS_delta_all_12 = [];
        CBS_theta_all_12 = [];
        CBS_alpha_all_12 = [];
        CBS_beta_all_12 = [];      
        CBS_delta_theta_all_12 = [];
        CBS_delta_alpha_all_12 = [];
        CBS_delta_beta_all_12 = [];
        CBS_theta_alpha_all_12 = [];
        CBS_theta_beta_all_12 = [];
        CBS_alpha_beta_all_12 = [];

        CBS_delta_all_21 = [];
        CBS_theta_all_21 = [];
        CBS_alpha_all_21 = [];
        CBS_beta_all_21 = [];
        CBS_delta_theta_all_21 = [];
        CBS_delta_alpha_all_21 = [];
        CBS_delta_beta_all_21 = [];
        CBS_theta_alpha_all_21 = [];
        CBS_theta_beta_all_21 = [];
        CBS_alpha_beta_all_21 = [];



%         mean_CBS_delta_all_21 = [];
%         mean_CBS_theta_all_21 = [];
%         mean_CBS_alpha_all_21 = [];
%         mean_CBS_beta_all_21 = [];
%         mean_CBS_beta1_all_21 = [];
%         mean_CBS_beta2_all_21 = [];
%         mean_CBS_beta3_all_21 = [];
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
            %compute Cross-bispectral 
            
            [Bspec1_ds,w] = BISPECDX(sig2_now, sig1_now, sig1_now, Fs, 6, 256);
            [Bspec2_ds,~] = BISPECDX(sig2_now, sig1_now, sig2_now, Fs, 6, 256);
            [freqs ,Bspec1] = bispect_sing(w,Bspec1_ds,Fs);
            [~ ,Bspec2] = bispect_sing(w,Bspec2_ds,Fs);

            delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
            theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
            alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
            beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
            beta1_idx = freqs >= beta_bd_1(1) & freqs <= beta_bd_1(2);
            beta2_idx = freqs >= beta_bd_2(1) & freqs <= beta_bd_2(2);
            beta3_idx = freqs >= beta_bd_3(1) & freqs <= beta_bd_3(2);

            CBS_delta_all_12 = [CBS_delta_all_12 mean(sum(Bspec1(delta_idx,delta_idx)))];
            CBS_theta_all_12 = [CBS_theta_all_12 mean(sum(Bspec1(theta_idx,theta_idx)))];
            CBS_alpha_all_12 = [CBS_alpha_all_12 mean(sum(Bspec1(alpha_idx,alpha_idx)))];
            CBS_beta_all_12 = [CBS_beta_all_12 mean(sum(Bspec1(beta_idx,beta_idx)))];
            CBS_delta_theta_all_12 = [CBS_delta_theta_all_12 mean(sum(Bspec1(delta_idx,theta_idx)))];
            CBS_delta_alpha_all_12 = [CBS_delta_alpha_all_12 mean(sum(Bspec1(delta_idx,alpha_idx)))];
            CBS_delta_beta_all_12 = [CBS_delta_beta_all_12 mean(sum(Bspec1(delta_idx,beta_idx)))];
            CBS_theta_alpha_all_12 = [CBS_theta_alpha_all_12 mean(sum(Bspec1(theta_idx,alpha_idx)))];
            CBS_theta_beta_all_12 = [CBS_theta_beta_all_12 mean(sum(Bspec1(theta_idx,beta_idx)))];
            CBS_alpha_beta_all_12 = [CBS_alpha_beta_all_12 mean(sum(Bspec1(alpha_idx,beta_idx)))];


            CBS_delta_all_21 = [CBS_delta_all_21 mean(sum(Bspec2(delta_idx,delta_idx)))];
            CBS_theta_all_21 = [CBS_theta_all_21 mean(sum(Bspec2(theta_idx,theta_idx)))];
            CBS_alpha_all_21 = [CBS_alpha_all_21 mean(sum(Bspec2(alpha_idx,alpha_idx)))];
            CBS_beta_all_21 = [CBS_beta_all_21 mean(sum(Bspec2(beta_idx,beta_idx)))];
            CBS_delta_theta_all_21 = [CBS_delta_theta_all_21 mean(sum(Bspec2(delta_idx,theta_idx)))];
            CBS_delta_alpha_all_21 = [CBS_delta_alpha_all_21 mean(sum(Bspec2(delta_idx,alpha_idx)))];
            CBS_delta_beta_all_21 = [CBS_delta_beta_all_21 mean(sum(Bspec2(delta_idx,beta_idx)))];
            CBS_theta_alpha_all_21 = [CBS_theta_alpha_all_21 mean(sum(Bspec2(theta_idx,alpha_idx)))];
            CBS_theta_beta_all_21 = [CBS_theta_beta_all_21 mean(sum(Bspec2(theta_idx,beta_idx)))];
            CBS_alpha_beta_all_21 = [CBS_alpha_beta_all_21 mean(sum(Bspec2(alpha_idx,beta_idx)))];
        end
        % add features into the results
        results = [results mean(CBS_delta_all_12)];
        results = [results mean(CBS_theta_all_12)];
        results = [results mean(CBS_alpha_all_12)];
        results = [results mean(CBS_beta_all_12)];
        results = [results mean(CBS_delta_theta_all_12)];
        results = [results mean(CBS_delta_alpha_all_12)];
        results = [results mean(CBS_delta_beta_all_12)];
        results = [results mean(CBS_theta_alpha_all_12)];
        results = [results mean(CBS_theta_beta_all_12)];
        results = [results mean(CBS_alpha_beta_all_12)];
        results = [results mean(CBS_delta_all_21)];
        results = [results mean(CBS_theta_all_21)];
        results = [results mean(CBS_alpha_all_21)];
        results = [results mean(CBS_beta_all_21)];
        results = [results mean(CBS_delta_theta_all_21)];
        results = [results mean(CBS_delta_alpha_all_21)];
        results = [results mean(CBS_delta_beta_all_21)];
        results = [results mean(CBS_theta_alpha_all_21)];
        results = [results mean(CBS_theta_beta_all_21)];
        results = [results mean(CBS_alpha_beta_all_21)];
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];
