function [results,dscrp,idx] = find_Bispectrum(y,Fs,ch_labels,split_length,results,dscrp)

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
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

for m = (1:ava_ch_num)

    sig1 = signal{m,1};
    sig1_ch = labels{m,1};

    % feature description
    dscrp{end+1} = string(sig1_ch)+'_BS_delta';
    dscrp{end+1} = string(sig1_ch)+'_BS_theta';
    dscrp{end+1} = string(sig1_ch)+'_BS_alpha';
    dscrp{end+1} = string(sig1_ch)+'_BS_beta';
    dscrp{end+1} = string(sig1_ch)+'_BS_delta_theta';
    dscrp{end+1} = string(sig1_ch)+'_BS_delta_alpha';
    dscrp{end+1} = string(sig1_ch)+'_BS_delta_beta';
    dscrp{end+1} = string(sig1_ch)+'_BS_theta_alpha';
    dscrp{end+1} = string(sig1_ch)+'_BS_theta_beta';
    dscrp{end+1} = string(sig1_ch)+'_BS_alpha_beta';

    % feature calculation container for each paired-channels
    BS_delta_all = [];
    BS_theta_all = [];
    BS_alpha_all = [];
    BS_beta_all = [];
    BS_delta_theta_all = [];
    BS_delta_alpha_all = [];
    BS_delta_beta_all = [];
    BS_theta_alpha_all = [];
    BS_theta_beta_all = [];
    BS_alpha_beta_all = [];


    % epoching the signal
    sig1_split = {};
    num_sig1 = floor(size(sig1,1)/split_len);

    % split the data into sub parts

    for i=0:1:num_sig1-1
        sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
    end

    sig1_split_size = size(sig1_split,1);
    min_sig_size = sig1_split_size;
    % traverse through the splited signals

    for k = 1:1:min_sig_size
        sig1_now = sig1_split{k};

        %feature calculation
        %compute Cross-bispectral

        [Bspec1_ds,w] = BISPECDX(sig1_now, sig1_now, sig1_now, Fs, 6, 256);
        [freqs ,Bspec1] = bispect_sing(w,Bspec1_ds,Fs);

        delta_idx = freqs >= delta_bd(1) & freqs <= delta_bd(2);
        theta_idx = freqs >= theta_bd(1) & freqs <= theta_bd(2);
        alpha_idx = freqs >= alpha_bd(1) & freqs <= alpha_bd(2);
        beta_idx = freqs >= beta_bd(1) & freqs <= beta_bd(2);
        beta1_idx = freqs >= beta_bd_1(1) & freqs <= beta_bd_1(2);
        beta2_idx = freqs >= beta_bd_2(1) & freqs <= beta_bd_2(2);
        beta3_idx = freqs >= beta_bd_3(1) & freqs <= beta_bd_3(2);

        BS_delta_all = [BS_delta_all mean(sum(Bspec1(delta_idx,delta_idx)))];
        BS_theta_all = [BS_theta_all mean(sum(Bspec1(theta_idx,theta_idx)))];
        BS_alpha_all = [BS_alpha_all mean(sum(Bspec1(alpha_idx,alpha_idx)))];
        BS_beta_all = [BS_beta_all mean(sum(Bspec1(beta_idx,beta_idx)))];
        BS_delta_theta_all = [BS_delta_theta_all mean(sum(Bspec1(delta_idx,theta_idx)))];
        BS_delta_alpha_all = [BS_delta_alpha_all mean(sum(Bspec1(delta_idx,alpha_idx)))];
        BS_delta_beta_all = [BS_delta_beta_all mean(sum(Bspec1(delta_idx,beta_idx)))];
        BS_theta_alpha_all = [BS_theta_alpha_all mean(sum(Bspec1(theta_idx,alpha_idx)))];
        BS_theta_beta_all = [BS_theta_beta_all mean(sum(Bspec1(theta_idx,beta_idx)))];
        BS_alpha_beta_all = [BS_alpha_beta_all mean(sum(Bspec1(alpha_idx,beta_idx)))];
    end
    % add features into the results
    results = [results mean(BS_delta_all)];
    results = [results mean(BS_theta_all)];
    results = [results mean(BS_alpha_all)];
    results = [results mean(BS_beta_all)];
    results = [results mean(BS_delta_theta_all)];
    results = [results mean(BS_delta_alpha_all)];
    results = [results mean(BS_delta_beta_all)];
    results = [results mean(BS_theta_alpha_all)];
    results = [results mean(BS_theta_beta_all)];
    results = [results mean(BS_alpha_beta_all)];
end
end_index = size(results,2);
idx = [start_index,end_index];
