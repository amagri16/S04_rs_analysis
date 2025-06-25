function [results,dscrp,idx] = find_Spearman(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/Junheng Li
split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
ava_ch_index = 1;
results = results;
dscrp = dscrp;

start_index = size(results,2)+1;
for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};


        % feature description
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Spearman';
        % feature calculation container for each paired-channels
        Spearman_all = [];

       
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1,1)/split_len);
        num_sig2 = floor(size(sig2,1)/split_len);
        % split the data into sub parts
        for i=0:1:num_sig1-1
            sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
        end

        for i=0:1:num_sig2-1
            sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
        end

        sig1_split_size = size(sig1_split,1);
        sig2_split_size = size(sig2_split,1);
        min_sig_size = min(sig1_split_size,sig2_split_size);
        % traverse through the splited signals

        for k = 1:1:min_sig_size
            % add feature analysis for each signal pearson
            sig1_now = sig1_split{k};
            sig2_now = sig2_split{k};
            Spearman_corr = corr(sig1_now, sig2_now, 'Type', 'Spearman');
            Spearman_all = [Spearman_all Spearman_corr];

        end
        results = [results mean(Spearman_all)]; 
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end

end_index = size(results,2);
idx = [start_index,end_index];