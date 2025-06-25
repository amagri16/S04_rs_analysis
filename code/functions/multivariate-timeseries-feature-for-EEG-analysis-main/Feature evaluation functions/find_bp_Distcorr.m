function [results,dscrp,idx] = find_bp_Distcorr(y,ch_labels,split_length,results,dscrp,filters)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
results = results;
dscrp = dscrp;
start_index = size(results,2)+1;
bands = ["delta","theta","alpha","beta1","beta2","beta3"];
%bands = ["delta","beta1"];

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};

        sig1_filt = [];
        sig2_filt = [];

        for i = 1:size(bands,2)
            sig1_filt(:,i) = filtfilt(filters{i}, sig1);
            sig2_filt(:,i) = filtfilt(filters{i}, sig2);
        end

        % epoching the signal
        sig1_split = {};
        sig2_split = {};
        num_sig1 = floor(size(sig1_filt,1)/split_len);
        num_sig2 = floor(size(sig2_filt,1)/split_len);
        


        % split the data into sub parts
        for i=0:1:num_sig1-1
            sig1_split{i+1,1}=sig1_filt(i*split_len+1:1:i*split_len+(split_len),:);
        end

        for i=0:1:num_sig2-1
            sig2_split{i+1,1}=sig2_filt(i*split_len+1:1:i*split_len+(split_len),:);
        end

        sig1_split_size = size(sig1_split,1);
        sig2_split_size = size(sig2_split,1);
        min_sig_size = min(sig1_split_size,sig2_split_size);

        for mm = 1:size(bands,2) % band-wise TE
            for nn = 1:size(bands,2)               

                dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_Distance_Corr'+'_'+bands(mm)+'_'+bands(nn);
                sig1_all = [];
                sig2_all = [];
                for i=0:1:num_sig1-1
                    sig1_all(:,i+1)=sig1_filt(i*split_len+1:1:i*split_len+(split_len),mm);
                end

                for i=0:1:num_sig2-1
                    sig2_all(:,i+1)=sig2_filt(i*split_len+1:1:i*split_len+(split_len),nn);
                end

               
                Corr_all = [];

                % traverse through the splited signals
                for k = 1:1:min_sig_size
                    sig1_now = sig1_split{k}(:,mm);
                    sig2_now = sig2_split{k}(:,nn);
                    
                    corre = distcorr(sig1_now, sig2_now);
                    Corr_all = [Corr_all corre];

                end
                results = [results mean(Corr_all)];
            end
        end
    end
end
end_index = size(results,2);
idx = [start_index,end_index];