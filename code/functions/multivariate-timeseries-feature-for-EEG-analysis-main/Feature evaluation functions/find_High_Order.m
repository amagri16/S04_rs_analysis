function [results,dscrp,idx] = find_High_Order(y,ch_labels,split_length,nplets,results,dscrp)
% Author: Tianyu Wei/Junheng Li
split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);

start_index = size(results,2)+1;

results = results;
dscrp = dscrp;
sig=[];
for i = 1:1:ava_ch_num
    sig = [sig signal{i}];
end

% feature description
dscrp{end+1} = 'mean_O_Information';
dscrp{end+1} = 'mean_S_Information';
dscrp{end+1} = 'mean_Redundancy';
dscrp{end+1} = 'mean_Synergy';
% feature calculation container for each epoch
O_all = [];
S_all = [];
Red_all = [];
Syn_all = [];
num_sig = floor(size(signal{1,1},1)/split_len);
for k = 1:1:num_sig
    % add feature analysis for each signal pearson
    sig_now = sig(1+split_len*(k-1):split_len*k,:)';
    [Oinfo,Sinfo,Red,Syn] = high_order(sig_now,nplets);
    O_all = [O_all mean(Oinfo)];
    S_all = [S_all mean(Sinfo)];
    Red_all = [Red_all Red];
    Syn_all = [Syn_all Syn];
end
results = [results mean(O_all)];
results = [results mean(S_all)];
results = [results mean(Red_all)];
results = [results mean(Syn_all)];
end_index = size(results,2);
idx = [start_index,end_index];
