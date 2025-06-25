function [results, dscrp, idx] = find_Spectral_Slope_avg(y,Fs,ch_labels,f_range,python_directory,iflinux,results,dscrp)

signal = y;
% how many channels exist
ava_ch_num = size(signal);
start_index = size(results,2)+1;

results_all = [];

for m = 1:ava_ch_num
    results_subj = [];
    y = signal{m,1}';
    y_ch_label = ch_labels{m,1};
    ft_dscrp_y = {};

    [psd_y,f] = pmtm(y,4,length(y),Fs);      %Multi-taper power-spectrum density analysis
    psd_y_norm = psd_y ./ sum(psd_y);

    exponent = extract_aperiodic_component(psd_y_norm, f, f_range, 0, python_directory, iflinux); % f_range should be in the format of [min_freq, max_freq]
    results_subj = [results_subj exponent];

    results_all = [results_all; results_subj];

    ft_dscrp_y = [ft_dscrp_y,{'Aperiodic_Component'}]; 
    
end

for i = 1:size(ft_dscrp_y,2)
    dscrp{end+1} = string(ft_dscrp_y{1,i});
end

results = [results mean(results_all,1)];

end_index = size(results,2);
idx = [start_index,end_index];
