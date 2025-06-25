function [results,dscrp,idx] = find_AMI(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/Junheng Li
split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal,1);
ava_ch_index = 1;
results = results;
dscrp = dscrp;
start_index = size(results,2)+1;
for m = (1:ava_ch_num)
%for m = ava_ch_num
    sig1 = signal{m,1};
    sig1_ch = labels{m,1};

    % feature description
    dscrp{end+1} = string(sig1_ch)+'_Automutual_Information';
    % feature calculation container for each paired-channels
    MI_all = [];

    sig1_split = {};
    num_sig1 = floor(size(sig1,1)/split_len);
    % split the data into sub parts
    for i=0:1:num_sig1-1
        sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
    end

    sig1_split_size = size(sig1_split,1);
    min_sig_size = sig1_split_size;
    % traverse through the splited signals

    for k = 1:1:min_sig_size-1
        sig1_now = sig1_split{k};
        sig2_now = sig1_split{k+1};
        MI = mutual_information(sig1_now', sig2_now', 'k', 4);
        MI_all = [MI_all mean(MI)];

    end
    results = [results mean(MI_all)];
    ava_ch_index = ava_ch_index+1;
end
end_index = size(results,2);
idx = [start_index,end_index];