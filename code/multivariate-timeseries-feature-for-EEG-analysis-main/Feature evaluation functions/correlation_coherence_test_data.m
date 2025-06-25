sig1 = sbj_eeg_combined.EEG_table.Signal{1,1};
sig2 = sbj_eeg_combined.EEG_table.Signal{2,1};
%% 
R = corrcoef(sig1,sig2);
noverlap = 256;
nfft = 512;
cxy = mscohere(sig1,sig2,hann(nfft),noverlap,nfft,256);
%% 
plot(cxy)
%% 
[R,P] = corrcoef(sig1(20000:30000),sig2(900000:910000))


%% split the signal into sub-parts
sig1_split = {};
sig2_split = {};
num_nfft_sig1 = floor(size(sig1,1)/nfft);
num_nfft_sig2 = floor(size(sig2,1)/nfft);
sig1_remain = [sig1(num_nfft_sig1*nfft+1:1:end,1)];
sig2_remain = [sig2(num_nfft_sig2*nfft+1:1:end,1)];


for i=0:1:num_nfft_sig1-1
    sig1_split{i+1,1}=sig1(i*nfft+1:1:i*nfft+(nfft));
end
if sig1_remain~=[]
sig1_split{num_nfft_sig1+1,1}=sig1_remain;
end

for i=0:1:num_nfft_sig2-1
    sig2_split{i+1,1}=sig2(i*nfft+1:1:i*nfft+(nfft));
end
if sig2_remain~=[]
sig2_split{num_nfft_sig2+1,1}=sig2_remain;
end

%% 
y = sig1(1:2048,1);
Fs = 256;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];

epc_len = 0.5;     % epoch length for feature extraction in 's'

len_y = round(length(y)/Fs);
num_epc_unit = len_y/epc_len;

for i = 1:num_epc_unit
    y_epc{i} = y((i-1)*Fs*epc_len+1:i*Fs*epc_len);
end
cnt = 0;
for i = 1:num_epc_unit-1
    for ii = i+1:num_epc_unit
        cnt= cnt+1;
        data1 = y_epc{i};
        data2 = y_epc{ii};
        [Cxy_y(cnt,:),f_cohere] = mscohere(data1,data2,[],[],0:0.5:30,Fs);
    end
end

%Cxy_y_avg = nanmean(Cxy_y);

Cxy_y_avg = mean(Cxy_y,'omitnan');

delta_idx = (f_cohere>delta_bd(1) & f_cohere<delta_bd(2));
theta_idx = (f_cohere>theta_bd(1) & f_cohere<theta_bd(2));
alpha_idx = (f_cohere>alpha_bd(1) & f_cohere<alpha_bd(2));
beta_idx = (f_cohere>beta_bd(1) & f_cohere<beta_bd(2));

% Temporal coherence in each frequency band
%ft.TCdelta = nanmean(Cxy_y_avg(delta_idx));
%ft.TCtheta = nanmean(Cxy_y_avg(theta_idx));
%ft.TCalpha = nanmean(Cxy_y_avg(alpha_idx));
%ft.TCbeta = nanmean(Cxy_y_avg(beta_idx));

ft.TCdelta = mean(Cxy_y_avg(delta_idx), 'omitnan');
ft.TCtheta = mean(Cxy_y_avg(theta_idx), 'omitnan');
ft.TCalpha = mean(Cxy_y_avg(alpha_idx), 'omitnan');
ft.TCbeta = mean(Cxy_y_avg(beta_idx), 'omitnan');



