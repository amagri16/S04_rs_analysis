function [results,dscrp] = find_DTF(y,Fs,results,dscrp)

signal = y;
Fs = 256;
signal = cell2mat(signal');
L = size(signal,1); % Number of samples
CH = size(signal,2); % Number of channels
EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};
N_freq = 600; % Number of frequency points
Fmax = 30;      % Maximum frequency limit in the PDC and DTF plots
win_len = Fs; % Window length
overlap = fix(win_len/2); % Overlap
start_point = 1;
s = 1;
Mode1 = 10; % ARFIT(10), Multichannel Yule-Walker(1) ---> BioSig toolbox 

[w, A_TI, C_TI, sbc, fpe, th] = arfit(signal, 1, 20, 'sbc'); % ---> ARFIT toolbox
[tmp,p_opt] = min(sbc);


while( start_point + win_len - 1 < L )
    seg = signal(start_point:start_point+win_len-1,:);
    seg = seg.*repmat(hamming(length(seg)),1,CH);
    [A_ST,RCF,C_ST] = mvar(seg, p_opt, Mode1); % ---> BioSig toolbox
    [PDC_ST(:,:,:,s), DTF_ST(:,:,:,s)] = PDC_DTF_matrix(A_ST,p_opt,Fs,Fmax,N_freq);
    start_point = start_point + (win_len-overlap);
    s = s + 1;
end


DTF_delta_1_2 = mean(mean(DTF_ST(1,2,1:floor(4*N_freq/30),:)));
DTF_delta_1_3 = mean(mean(DTF_ST(1,3,1:floor(4*N_freq/30),:)));
DTF_delta_2_3 = mean(mean(DTF_ST(2,3,1:floor(4*N_freq/30),:)));

DTF_theta_1_2 = mean(mean(DTF_ST(1,2,floor(4*N_freq/30):floor(8*N_freq/30),:)));
DTF_theta_1_3 = mean(mean(DTF_ST(1,3,floor(4*N_freq/30):floor(8*N_freq/30),:)));
DTF_theta_2_3 = mean(mean(DTF_ST(2,3,floor(4*N_freq/30):floor(8*N_freq/30),:)));

DTF_alpha_1_2 = mean(mean(DTF_ST(1,2,floor(8*N_freq/30):floor(12*N_freq/30),:)));
DTF_alpha_1_3 = mean(mean(DTF_ST(1,3,floor(8*N_freq/30):floor(12*N_freq/30),:)));
DTF_alpha_2_3 = mean(mean(DTF_ST(2,3,floor(8*N_freq/30):floor(12*N_freq/30),:)));

DTF_sigma_1_2 = mean(mean(DTF_ST(1,2,floor(12*N_freq/30):floor(18*N_freq/30),:)));
DTF_sigma_1_3 = mean(mean(DTF_ST(1,3,floor(12*N_freq/30):floor(18*N_freq/30),:)));
DTF_sigma_2_3 = mean(mean(DTF_ST(2,3,floor(12*N_freq/30):floor(18*N_freq/30),:)));

dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_DTF_delta';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_DTF_delta';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_DTF_delta';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_DTF_theta';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_DTF_theta';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_DTF_theta';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_DTF_alpha';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_DTF_alpha';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_DTF_alpha';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_DTF_sigma';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_DTF_sigma';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_DTF_sigma';
results = [results DTF_delta_1_2];
results = [results DTF_delta_1_3];
results = [results DTF_delta_2_3];
results = [results DTF_theta_1_2];
results = [results DTF_theta_1_3];
results = [results DTF_theta_2_3];
results = [results DTF_alpha_1_2];
results = [results DTF_alpha_1_3];
results = [results DTF_alpha_2_3];
results = [results DTF_sigma_1_2];
results = [results DTF_sigma_1_3];
results = [results DTF_sigma_2_3];

%%
PDC_delta_1_2 = mean(mean(PDC_ST(1,2,1:floor(4*N_freq/30),:)));
PDC_delta_1_3 = mean(mean(PDC_ST(1,3,1:floor(4*N_freq/30),:)));
PDC_delta_2_3 = mean(mean(PDC_ST(2,3,1:floor(4*N_freq/30),:)));

PDC_theta_1_2 = mean(mean(PDC_ST(1,2,floor(4*N_freq/30):floor(8*N_freq/30),:)));
PDC_theta_1_3 = mean(mean(PDC_ST(1,3,floor(4*N_freq/30):floor(8*N_freq/30),:)));
PDC_theta_2_3 = mean(mean(PDC_ST(2,3,floor(4*N_freq/30):floor(8*N_freq/30),:)));

PDC_alpha_1_2 = mean(mean(PDC_ST(1,2,floor(8*N_freq/30):floor(12*N_freq/30),:)));
PDC_alpha_1_3 = mean(mean(PDC_ST(1,3,floor(8*N_freq/30):floor(12*N_freq/30),:)));
PDC_alpha_2_3 = mean(mean(PDC_ST(2,3,floor(8*N_freq/30):floor(12*N_freq/30),:)));

PDC_sigma_1_2 = mean(mean(PDC_ST(1,2,floor(12*N_freq/30):floor(18*N_freq/30),:)));
PDC_sigma_1_3 = mean(mean(PDC_ST(1,3,floor(12*N_freq/30):floor(18*N_freq/30),:)));
PDC_sigma_2_3 = mean(mean(PDC_ST(2,3,floor(12*N_freq/30):floor(18*N_freq/30),:)));

dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_PDC_delta';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_PDC_delta';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_PDC_delta';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_PDC_theta';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_PDC_theta';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_PDC_theta';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_PDC_alpha';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_PDC_alpha';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_PDC_alpha';
dscrp{end+1} = 'Fz-Cz_to_Cz-Oz_PDC_sigma';
dscrp{end+1} = 'Fz-Cz_to_C4-M1_PDC_sigma';
dscrp{end+1} = 'Cz-Oz_to_C4-M1_PDC_sigma';
results = [results PDC_delta_1_2];
results = [results PDC_delta_1_3];
results = [results PDC_delta_2_3];
results = [results PDC_theta_1_2];
results = [results PDC_theta_1_3];
results = [results PDC_theta_2_3];
results = [results PDC_alpha_1_2];
results = [results PDC_alpha_1_3];
results = [results PDC_alpha_2_3];
results = [results PDC_sigma_1_2];
results = [results PDC_sigma_1_3];
results = [results PDC_sigma_2_3];

%% 
DTF_delta_2_1 = mean(mean(DTF_ST(2,1,1:floor(4*N_freq/30),:)));
DTF_delta_3_1 = mean(mean(DTF_ST(3,1,1:floor(4*N_freq/30),:)));
DTF_delta_3_2 = mean(mean(DTF_ST(3,2,1:floor(4*N_freq/30),:)));

DTF_theta_2_1 = mean(mean(DTF_ST(2,1,floor(4*N_freq/30):floor(8*N_freq/30),:)));
DTF_theta_3_1 = mean(mean(DTF_ST(3,1,floor(4*N_freq/30):floor(8*N_freq/30),:)));
DTF_theta_3_2 = mean(mean(DTF_ST(3,2,floor(4*N_freq/30):floor(8*N_freq/30),:)));

DTF_alpha_2_1 = mean(mean(DTF_ST(2,1,floor(8*N_freq/30):floor(12*N_freq/30),:)));
DTF_alpha_3_1 = mean(mean(DTF_ST(3,1,floor(8*N_freq/30):floor(12*N_freq/30),:)));
DTF_alpha_3_2 = mean(mean(DTF_ST(3,2,floor(8*N_freq/30):floor(12*N_freq/30),:)));

DTF_sigma_2_1 = mean(mean(DTF_ST(2,1,floor(12*N_freq/30):floor(18*N_freq/30),:)));
DTF_sigma_3_1 = mean(mean(DTF_ST(3,1,floor(12*N_freq/30):floor(18*N_freq/30),:)));
DTF_sigma_3_2 = mean(mean(DTF_ST(3,2,floor(12*N_freq/30):floor(18*N_freq/30),:)));

dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_DTF_delta';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_DTF_delta';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_DTF_delta';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_DTF_theta';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_DTF_theta';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_DTF_theta';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_DTF_alpha';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_DTF_alpha';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_DTF_alpha';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_DTF_sigma';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_DTF_sigma';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_DTF_sigma';

results = [results DTF_delta_2_1];
results = [results DTF_delta_3_1];
results = [results DTF_delta_3_2];
results = [results DTF_theta_2_1];
results = [results DTF_theta_3_1];
results = [results DTF_theta_3_2];
results = [results DTF_alpha_2_1];
results = [results DTF_alpha_3_1];
results = [results DTF_alpha_3_2];
results = [results DTF_sigma_2_1];
results = [results DTF_sigma_3_1];
results = [results DTF_sigma_3_2];

%%
PDC_delta_2_1 = mean(mean(PDC_ST(2,1,1:floor(4*N_freq/30),:)));
PDC_delta_3_1 = mean(mean(PDC_ST(3,1,1:floor(4*N_freq/30),:)));
PDC_delta_3_2 = mean(mean(PDC_ST(3,2,1:floor(4*N_freq/30),:)));

PDC_theta_2_1 = mean(mean(PDC_ST(2,1,floor(4*N_freq/30):floor(8*N_freq/30),:)));
PDC_theta_3_1 = mean(mean(PDC_ST(3,1,floor(4*N_freq/30):floor(8*N_freq/30),:)));
PDC_theta_3_2 = mean(mean(PDC_ST(3,2,floor(4*N_freq/30):floor(8*N_freq/30),:)));

PDC_alpha_2_1 = mean(mean(PDC_ST(2,1,floor(8*N_freq/30):floor(12*N_freq/30),:)));
PDC_alpha_3_1 = mean(mean(PDC_ST(3,1,floor(8*N_freq/30):floor(12*N_freq/30),:)));
PDC_alpha_3_2 = mean(mean(PDC_ST(3,2,floor(8*N_freq/30):floor(12*N_freq/30),:)));

PDC_sigma_2_1 = mean(mean(PDC_ST(2,1,floor(12*N_freq/30):floor(18*N_freq/30),:)));
PDC_sigma_3_1 = mean(mean(PDC_ST(3,1,floor(12*N_freq/30):floor(18*N_freq/30),:)));
PDC_sigma_3_2 = mean(mean(PDC_ST(3,2,floor(12*N_freq/30):floor(18*N_freq/30),:)));

dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_PDC_delta';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_PDC_delta';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_PDC_delta';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_PDC_theta';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_PDC_theta';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_PDC_theta';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_PDC_alpha';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_PDC_alpha';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_PDC_alpha';
dscrp{end+1} = 'Cz-Oz_to_Fz-Cz_PDC_sigma';
dscrp{end+1} = 'C4-M1_to_Fz-Cz_PDC_sigma';
dscrp{end+1} = 'C4-M1_to_Cz-Oz_PDC_sigma';

results = [results PDC_delta_2_1];
results = [results PDC_delta_3_1];
results = [results PDC_delta_3_2];
results = [results PDC_theta_2_1];
results = [results PDC_theta_3_1];
results = [results PDC_theta_3_2];
results = [results PDC_alpha_2_1];
results = [results PDC_alpha_3_1];
results = [results PDC_alpha_3_2];
results = [results PDC_sigma_2_1];
results = [results PDC_sigma_3_1];
results = [results PDC_sigma_3_2];
