%% correlation and coherence of the EEG signal
%function correlation_coherence, c_value = correlation_coherence(patientObj,EEG_name,f_s);

% mscohere
% Magnitude-squared coherence mscohere estimates the magnitude-squared
% coherence function using Welch s overlapped averaged periodogram method 

% corrcoef returns the matrix of correlation coefficients

% binary classification
% https://uk.mathworks.com/help/matlab/ref/xcorr.html
% build a struct
% corre 1/2 peak value length of shift
% coherence -4 frequencies -average
% use xcov
%% 
%clear,clc;
load('eeg_combined_example.mat')
patientObj = sbj_eeg_combined;
EEG_name = 'Onset_period';

%% Check if subject is discarded

if patientObj.Sbj_discard
    disp('Subject discarded, no action will be done');
    return
end

% Check if there are channels available to use
if sum(1-[patientObj.EEG_table.discarded]) == 0
    patientObj.Sbj_discard = 1;
    patientObj.Discard_info = 'All EEG channels discarded';
    disp(['All EEG channels dicarded, no processing for subject No.',...
        num2str(patientObj.patient_original_id)])
    return
end

%% find all channels that are available

EEG_table_height = height(patientObj.EEG_table);
ava_ch = [];

for i = 1:1:EEG_table_height
    if ~patientObj.EEG_table.discarded(i)
        ava_ch = [ava_ch i];
    end
end

%% Cross-compute results for all pairs of signal

num_ava_ch = size(ava_ch);

% the ava_ch is not always continous it might be [1 2 4 5 6 7 8 10]
% using both index and value corresponds with the index
f_s = 256;
noverlap = 256;
split_len = 1024;
ava_ch_index = 1;
nfft = 0:0.5:30;

delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];


%% split the signal into sub-parts

% this block of code evenly divide the signal into a number of epoches with
% size equal to the number of data points taken when computing DFT(the DFT
% of the autocorrelation function when finding the coherence), then proceed
% with computing the correlation between each epoch and compile the results
% into a cell. This could makes the plot easier.

ava_ch_index = 1;
corr_results = {};
cxy_results = {};
for m = ava_ch(1:end-1)
    for n = ava_ch(ava_ch_index+1:1:end)
        sig1 = patientObj.EEG_table.(EEG_name){m,1};
        sig2 = patientObj.EEG_table.(EEG_name){n,1};
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1,1)/split_len);
        num_sig2 = floor(size(sig2,1)/split_len);
        sig1_remain = [sig1(num_sig1*split_len+1:1:end,1)];
        sig2_remain = [sig2(num_sig2*split_len+1:1:end,1)];

        for i=0:1:num_sig1-1
            sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
        end

        for i=0:1:num_sig2-1
            sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
        end
        
        corre_i_j = {};
        cxy_i_j = {};
        corre_i_j_index = 1;
        cxy_i_j_index = 1;
        sig1_split_size = size(sig1_split,1);
        sig2_split_size = size(sig2_split,1);
        min_sig_size = min(sig1_split_size,sig2_split_size);
        
        % traverse through the splited signals
        for k = 1:1:min_sig_size
            R = xcorr(sig1_split{k,1},sig2_split{k,1});
            corre_i_j{k,1} = mean(R);
            corre_i_j_index = corre_i_j_index+1;
            cxy_i_j{k,1} = mscohere(sig1_split{k,1},sig2_split{k,1},[],[],nfft,f_s);
            cxy_i_j_index = cxy_i_j_index+1;
        end

        corr_results{m,n} = mean(transpose(cell2mat(corre_i_j)));
        cxy_results{m,n} = mean(transpose(cell2mat(cxy_i_j)),2);
        corr_dscrp{m,n} = 'cor_'+string(m)+'_'+string(n);
        cxy_dscrp{m,n} = 'coh_'+string(m)+'_'+string(n);
    end
    ava_ch_index = ava_ch_index + 1;
end
%% coherence on each band
for i=1:1:m
    for j=1:1:n
        if isempty(cxy_results{i,j})
            cxy_results{i,j} = nan;
        else
            delta_idx = (nfft>delta_bd(1) & nfft<delta_bd(2));
            theta_idx = (nfft>theta_bd(1) & nfft<theta_bd(2));
            alpha_idx = (nfft>alpha_bd(1) & nfft<alpha_bd(2));
            beta_idx = (nfft>beta_bd(1) & nfft<beta_bd(2));

            mean_delta_coh = mean(cxy_results{i,j}(delta_idx));
            mean_theta_coh = mean(cxy_results{i,j}(theta_idx));
            mean_alpha_coh = mean(cxy_results{i,j}(alpha_idx));
            mean_beta_coh = mean(cxy_results{i,j}(beta_idx));
            cxy_results{i,j} = [mean_delta_coh, mean_theta_coh, mean_alpha_coh, mean_beta_coh];
        end
    end
end

%% add the tables into the patientObj
patientObj.EEG_table.('Corr_Results_'+string(EEG_name)){1} = corr_results;
patientObj.EEG_table.('Coherence_Results_'+string(EEG_name)){1} = cxy_results;
%% coherence test
test_sig1 = sig1(7500:8000,1);
test_sig2 = sig2(7500:8000,1);
[cxy_test,f_cohere] = mscohere(test_sig1,test_sig2,[],[],nfft,256);
[coh_max,coh_max_index] = max(cxy_test);
f_cohere_max = f_cohere(coh_max_index);

delta_idx = (f_cohere>delta_bd(1) & f_cohere<delta_bd(2));
theta_idx = (f_cohere>theta_bd(1) & f_cohere<theta_bd(2));
alpha_idx = (f_cohere>alpha_bd(1) & f_cohere<alpha_bd(2));
beta_idx = (f_cohere>beta_bd(1) & f_cohere<beta_bd(2));

mean_delta_coh = mean(cxy_test(delta_idx), 'omitnan');
mean_theta_coh = mean(cxy_test(theta_idx), 'omitnan');
mean_alpha_coh = mean(cxy_test(alpha_idx), 'omitnan');
mean_beta_coh = mean(cxy_test(beta_idx), 'omitnan');

%% correlation test
corr_test = xcorr(test_sig1,test_sig2);
[corr_max,corr_max_index] = max(corr_test);
half_max = 0.5*corr_max;
corr_idx = (corr_test>half_max);

[labeledX, numRegions] = bwlabel(corr_idx);
props = regionprops(labeledX, 'Area', 'PixelList');
regionLengths = [props.Area];
half_lag = find(labeledX==labeledX(corr_max_index));
half_lag_len = size(tt);


