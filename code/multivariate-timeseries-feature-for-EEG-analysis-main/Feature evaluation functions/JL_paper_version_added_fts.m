function [results, dscrp, idx] = JL_paper_version_added_fts(y,Fs,ch_labels,f_range,python_directory,iflinux,results,dscrp)

signal = y;
% how many channels exist
ava_ch_num = size(signal);
start_index = size(results,2)+1;

sigma_bd = [12,16];
   

for m = 1:ava_ch_num
    y = signal{m,1}';

    y_MA = y;
    y = y_MA(2:end) - y_MA(1:end-1);
    y = [y,0];

    order = 2;
    lowFreq0 = 0.1;    % 0.1
    highFreq0 = 32;
    [b,a] = butter(order, [lowFreq0 highFreq0]/(Fs/2), 'bandpass');
    y = filtfilt(b,a,y);

    y_ch_label = ch_labels{m,1};
    ft_dscrp_y = {};

    [psd_y,f] = pmtm(y,4,length(y),Fs);      %Multi-taper power-spectrum density analysis
    psd_y_norm = psd_y ./ sum(psd_y);    % Normalize to total power
    sigma_idx = (f>=sigma_bd(1) & f<sigma_bd(2));

    sigma_power = mean(psd_y_norm(sigma_idx),'omitmissing');


    results = [results sigma_power];
    ft_dscrp_y = [ft_dscrp_y,{'sigma'}]; 

    for i = 1:size(ft_dscrp_y,2)
        dscrp{end+1} = string(y_ch_label)+'_'+string(ft_dscrp_y{1,i});       
    end
end


for m = 1:ava_ch_num
    y = signal{m,1}';
    y_ch_label = ch_labels{m,1};
    ft_dscrp_y = {};


    binarized_y = y > mean(y,2);
    lempel_ziv = LZ76(binarized_y);
    entropy_rate_lz = lempel_ziv*log2(length(binarized_y))/length(binarized_y);

    results = [results entropy_rate_lz];

    ft_dscrp_y = [ft_dscrp_y,{'LZComplexity'}];

    for i = 1:size(ft_dscrp_y,2)
        dscrp{end+1} = string(y_ch_label)+'_'+string(ft_dscrp_y{1,i});       
    end
end



for m = 1:ava_ch_num
    y = signal{m,1}';

    y_ch_label = ch_labels{m,1};
    ft_dscrp_y = {};

    [psd_y,f] = pmtm(y,4,length(y),Fs);      %Multi-taper power-spectrum density analysis
    psd_y_norm = psd_y ./ sum(psd_y);

    exponent = extract_aperiodic_component(psd_y_norm, f, f_range, 0, python_directory, iflinux); % f_range should be in the format of [min_freq, max_freq]
    results = [results exponent];
    ft_dscrp_y = [ft_dscrp_y,{'Exp1f'}]; 

    for i = 1:size(ft_dscrp_y,2)
        dscrp{end+1} = string(y_ch_label)+'_'+string(ft_dscrp_y{1,i});       
    end
end

end_index = size(results,2);
idx = [start_index,end_index];
