function [results, dscrp, idx] = find_Complexity(y,ch_labels,results,dscrp)

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

    entropy_rate = CTWEntropyRate(y);
    binarized_y = y > mean(y,2);
    lempel_ziv = LZ76(binarized_y);
    entropy_rate_lz = lempel_ziv*log2(length(binarized_y))/length(binarized_y);

    results_subj = [results_subj, entropy_rate];
    results_subj = [results_subj, lempel_ziv];
    results_subj = [results_subj, entropy_rate_lz];
    
    results_all = [results_all; results_subj];

    ft_dscrp_y = [ft_dscrp_y,{'Entropy_Rate_CTW','Lempel_Ziv', 'Entropy_Rate_LZ'}];
end

for i = 1:size(ft_dscrp_y,2)
    dscrp{end+1} = string(ft_dscrp_y{1,i});
end

results = [results mean(results_all,1)];

end_index = size(results,2);
idx = [start_index,end_index];
