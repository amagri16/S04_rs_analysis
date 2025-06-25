function [results,dscrp,idx] = find_TE_corrected(y,ch_labels,split_length,results,dscrp)
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
for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Transfer_Entropy_PTE';
        dscrp{end+1} = string(sig2_ch)+'_vs_'+string(sig1_ch)+'_Transfer_Entropy_PTE';
        % feature calculation container for each paired-channels
        TE1_all = [];
        TE2_all = [];
        % epoching the signal
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
            sig1_now = sig1_split{k};
            sig2_now = sig2_split{k};
            % Compute the Spectral Entropy
            [~,eLag,eDim] = phaseSpaceReconstruction([sig1_now,sig2_now]);
            TE_result = TE(sig1_now,sig2_now,4,0,eDim,eLag,0);
            TE1_all = [TE1_all; TE_result(1,1)];
            TE2_all = [TE2_all; TE_result(1,2)];

        end
        results = [results mean(TE1_all)];
        results = [results mean(TE2_all)];
        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];