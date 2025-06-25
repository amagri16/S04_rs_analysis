function [results, dscrp, idx] = find_DTF(y,Fs,channel_labels,results,dscrp)

signal = y;
signal = cell2mat(signal');

% signal = signal(2:end,:) - signal(1:end-1,:); % Apply differentiation filter

L = size(signal,1); % Number of samples
CH = size(signal,2); % Number of channels
EEG_chs = channel_labels;
N_freq = 256; % Number of frequency points
Fmax = 30;      % Maximum frequency limit in the PDC and DTF plots
win_len = Fs; % Window length
overlap = fix(win_len/2); % Overlap
start_point = 1;
s = 1;
Mode1 = 10; % ARFIT(10), Multichannel Yule-Walker(1) ---> BioSig toolbox 

[w, A_TI, C_TI, sbc, fpe, th] = arfit(signal, 1, 20, 'sbc'); % ---> ARFIT toolbox
[tmp,p_opt] = min(sbc);
start_index = size(results,2)+1;

while( start_point + win_len - 1 < L )
    seg = signal(start_point:start_point+win_len-1,:);
    seg = seg.*repmat(hamming(length(seg)),1,CH);
    [A_ST,RCF,C_ST] = mvar(seg, p_opt, Mode1); % ---> BioSig toolbox
    [PDC_ST(:,:,:,s), DTF_ST(:,:,:,s)] = PDC_DTF_matrix(A_ST,p_opt,Fs,Fmax,N_freq);
    start_point = start_point + (win_len-overlap);
    s = s + 1;
end

% Define frequency bands and labels
freq_bands = [1, 4; 4, 8; 8, 12; 12, 18; 18, 24; 24, 30];
band_labels = {'delta', 'theta', 'alpha', 'beta1','beta2','beta3'};

% Compute the mean values, update the description and results arrays
for i = 1:CH
    for j = 1:CH
        if i ~= j
            for k = 1:size(freq_bands, 1)
                f_start = floor(freq_bands(k, 1) * N_freq / Fmax);
                f_end = floor(freq_bands(k, 2) * N_freq / Fmax);
                
                DTF_val = mean(mean(DTF_ST(i, j, f_start:f_end, :)));
                PDC_val = mean(mean(PDC_ST(i, j, f_start:f_end, :)));
                
                dscrp{end+1} = [EEG_chs{i}, '_to_', EEG_chs{j}, '_DTF_', band_labels{k}];
                dscrp{end+1} = [EEG_chs{i}, '_to_', EEG_chs{j}, '_PDC_', band_labels{k}];
                results = [results DTF_val PDC_val];
            end
        end
    end
end
end_index = size(results,2);
idx = [start_index,end_index];

