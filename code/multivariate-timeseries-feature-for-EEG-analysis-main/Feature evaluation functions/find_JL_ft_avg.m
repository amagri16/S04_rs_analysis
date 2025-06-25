function [results, dscrp, idx] = find_JL_ft_avg(y,Fs,ch_labels,results,dscrp,filters)

signal = y;
% how many channels exist
ava_ch_num = size(signal);
start_index = size(results,2)+1;

results_all = [];
for m = 1:ava_ch_num
    results_subj = [];

    % Band frequency specification
    y = signal{m,1};
    y_ch_label = ch_labels{m,1};

    if size(y,2)<size(y,1)     %make sure it's a row vector
        y = transpose(y);
    end


    %% PSD estimation part

    delta_bd = [1,4];
    theta_bd = [4,8];
    alpha_bd = [8,12];
    beta_bd = [13,30];
    beta_narrow_bd = [25,29];
    gamma_bd = [30,35];

     ft_dscrp_y = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency band features
    ifproc = 1;

    % y_MA = y;
    % y = y_MA(2:end) - y_MA(1:end-1); 
    % y = [y,0];

    if ifproc     % Normal pre-processing
        % filter EEG, broad BP
        order = 2;
        lowFreq0 = 0.1;    % 0.1
        highFreq0 = 32;
        [b,a] = butter(order, [lowFreq0 highFreq0]/(Fs/2), 'bandpass'); %pre-filtering via Butterworth bandpass
        y_proc = filtfilt(b,a,y);

    else
        y_proc = y;
    end


    [psd_y,f] = pmtm(y_proc,4,length(y_proc),Fs);      %Multi-taper power-spectrum density analysis

    psd_y_norm = psd_y ./ sum(psd_y);    % Normalize to total power

    delta_idx = (f>delta_bd(1) & f<delta_bd(2));
    theta_idx = (f>theta_bd(1) & f<theta_bd(2));
    alpha_idx = (f>alpha_bd(1) & f<alpha_bd(2));
    beta_idx = (f>beta_bd(1) & f<beta_bd(2));
    beta_narrow_idx = (f>beta_narrow_bd(1) & f<beta_narrow_bd(2));
    gamma_idx = (f>gamma_bd(1) & f<gamma_bd(2));

    delta_power = mean(psd_y_norm(delta_idx),'omitnan');
    theta_power = mean(psd_y_norm(theta_idx),'omitnan');
    alpha_power = mean(psd_y_norm(alpha_idx),'omitnan');
    beta_power = mean(psd_y_norm(beta_idx),'omitnan');
    beta_narrow_power = mean(psd_y_norm(beta_narrow_idx),'omitnan');
    gamma_power = mean(psd_y_norm(gamma_idx),'omitnan');

    results_subj = [results_subj log(delta_power)];
    results_subj = [results_subj log(theta_power)];
    results_subj = [results_subj log(alpha_power)];
    results_subj = [results_subj log(beta_power)];
    results_subj = [results_subj log(beta_narrow_power)];
    results_subj = [results_subj log(gamma_power)];
    results_subj = [results_subj log(delta_power/theta_power)];
    results_subj = [results_subj log(delta_power/alpha_power)];
    results_subj = [results_subj log(delta_power/beta_power)];
    results_subj = [results_subj log(theta_power/alpha_power)];
    results_subj = [results_subj log(theta_power/beta_power)];
    results_subj = [results_subj log(alpha_power/beta_power)];

    results_subj = [results_subj log((delta_power+theta_power)/(alpha_power+beta_power))];    % Transition power ratio

    ft_dscrp_y = [ft_dscrp_y,{'delta','theta','alpha','beta 25-29','beta','gamma','DTratio','DAratio','DBratio','TAratio','TBratio','ABratio','Transratio'}];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Peak Oscillation frequency
    f_delta = f(delta_idx);
    results_subj = [results_subj f_delta(find(psd_y_norm(delta_idx) == max(psd_y_norm(delta_idx)),1))];

    f_theta = f(theta_idx);
    results_subj = [results_subj f_theta(find(psd_y_norm(theta_idx) == max(psd_y_norm(theta_idx)),1))];

    f_alpha = f(alpha_idx);
    results_subj = [results_subj f_alpha(find(psd_y_norm(alpha_idx) == max(psd_y_norm(alpha_idx)),1))];

    f_beta = f(beta_idx);
    results_subj = [results_subj f_beta(find(psd_y_norm(beta_idx) == max(psd_y_norm(beta_idx)),1))];

    f_beta_narrow = f(beta_narrow_idx);
    results_subj = [results_subj f_beta_narrow(find(psd_y_norm(beta_narrow_idx) == max(psd_y_norm(beta_narrow_idx)),1))];

    f_gamma = f(gamma_idx);
    results_subj = [results_subj f_gamma(find(psd_y_norm(gamma_idx) == max(psd_y_norm(gamma_idx)),1))];

    ft_dscrp_y = [ft_dscrp_y,{'deltapeak','thetapeak','alphapeak','betapeak','beta 25-29 peak','gammapeak'}];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temporal coherence features

    epc_len = 0.5;     % epoch length for feature extraction in 's'

    len_y = round(length(y)/Fs);
    num_epc_unit = len_y/epc_len-1;

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

    Cxy_y_avg = mean(Cxy_y,'omitnan');

    delta_idx = (f_cohere>delta_bd(1) & f_cohere<delta_bd(2));
    theta_idx = (f_cohere>theta_bd(1) & f_cohere<theta_bd(2));
    alpha_idx = (f_cohere>alpha_bd(1) & f_cohere<alpha_bd(2));
    beta_idx = (f_cohere>beta_bd(1) & f_cohere<beta_bd(2));

    results_subj = [results_subj mean(Cxy_y_avg(delta_idx), 'omitnan')];
    results_subj = [results_subj mean(Cxy_y_avg(theta_idx), 'omitnan')];
    results_subj = [results_subj mean(Cxy_y_avg(alpha_idx), 'omitnan')];
    results_subj = [results_subj mean(Cxy_y_avg(beta_idx), 'omitnan')];

    ft_dscrp_y = [ft_dscrp_y,{'TCdelta','TCtheta','TCalpha','TCbeta'}];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Phase-amplitude coupling features
    % Filter chosen is zero-phase, even order FIR filter;
    y_filt_delta = filtfilt(filters{1},y);
    y_filt_theta = filtfilt(filters{2},y);
    y_filt_alpha = filtfilt(filters{3},y);
    y_filt_beta = filtfilt(filters{4},y);

    y_delta_analytic = hilbert(y_filt_delta);
    y_theta_analytic = hilbert(y_filt_theta);
    y_alpha_analytic = hilbert(y_filt_alpha);
    y_beta_analytic = hilbert(y_filt_beta);

    % Amplitudes
    amp_theta = abs(y_theta_analytic);
    amp_alpha = abs(y_alpha_analytic);
    amp_beta = abs(y_beta_analytic);

    % Phases
    phase_delta = angle(y_delta_analytic);
    phase_theta = angle(y_theta_analytic);
    phase_alpha = angle(y_alpha_analytic);

    idx_cycle_delta = find_cycles(phase_delta);
    idx_cycle_theta = find_cycles(phase_theta);
    idx_cycle_alpha = find_cycles(phase_alpha);

    zt_DT = amp_theta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
    zt_DA = amp_alpha(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
    zt_DB = amp_beta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
    zt_TA = amp_alpha(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
    zt_TB = amp_beta(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
    zt_AB = amp_beta(1:idx_cycle_alpha).*exp(j*phase_alpha(1:idx_cycle_alpha));

    results_subj = [results_subj abs(sum(zt_DT)/sum(amp_theta(1:idx_cycle_delta)))];
    results_subj = [results_subj abs(sum(zt_DA)/sum(amp_alpha(1:idx_cycle_delta)))];
    results_subj = [results_subj abs(sum(zt_DB)/sum(amp_beta(1:idx_cycle_delta)))];
    results_subj = [results_subj abs(sum(zt_TA)/sum(amp_alpha(1:idx_cycle_theta)))];
    results_subj = [results_subj abs(sum(zt_TB)/sum(amp_beta(1:idx_cycle_theta)))];
    results_subj = [results_subj abs(sum(zt_AB)/sum(amp_beta(1:idx_cycle_alpha)))];

    ft_dscrp_y = [ft_dscrp_y,{'DTcouple','DAcouple','DBcouple','TAcouple','TBcouple','ABcouple'}];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Catch-22 cannonical features

    if size(y,2)>size(y,1)     %make sure it's a column vector
        y = transpose(y);
    end

    [ft_values,ft_name] = catch22_all(y); % run catch22-master/wrap_Matlab/mexAll.m before executing this line
    for i = 1:length(ft_values)
        results_subj = [results_subj ft_values(i)];
    end
    ft_dscrp_y = [ft_dscrp_y,ft_name];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    results_all = [results_all; results_subj];

end

for i = 1:size(ft_dscrp_y,2)
    dscrp{end+1} = string(ft_dscrp_y{1,i});
end

results = [results mean(results_all,1)];

end_index = size(results,2);
idx = [start_index,end_index];

function idx = find_cycles(phase)
stp = phase(1);
for i = length(phase):-1:2
    if phase(2)>stp
        if phase(i) >= stp && phase(i-1)<= stp
            idx = i;
            return
        end
    end
    if phase(2)<stp
        if phase(i) <= stp && phase(i-1)>= stp
            idx = i;
            return
        end
    end
end