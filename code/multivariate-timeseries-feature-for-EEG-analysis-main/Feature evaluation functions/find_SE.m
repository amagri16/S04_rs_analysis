function [results,dscrp,idx] = find_SE(y,Fs,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
results = results;
dscrp = dscrp;
start_index = size(results,2)+1;
ts = 1/256;
for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spectral_Entropy';

        % Compute the Spectral Entropy
        SE1 = pentropy(sig1,ts);
        SE2 = pentropy(sig2,ts);

        % add features into the results
        SE = corr(SE1,SE2,'Type', 'Spearman');

        results = [results SE];

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];