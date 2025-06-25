function [results,dscrp,idx] = find_TE_final(y,ch_labels,split_length,results,dscrp)
% Author: Tianyu Wei/ Junheng Li

split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal,1);
ava_ch_index = 1;
results = results;
dscrp = dscrp;
start_index = size(results,2)+1;
lagSet = -2 : 2; %WAS -5 TO 5
for m = (1:ava_ch_num-1)
    % for n = (m+1:1:ava_ch_num)
    for n = ava_ch_num
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        
        % feature description
        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Transfer_Entropy';
        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Transfer_Entropy';
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


        sig1_all = [];
        sig2_all = [];
        for i=0:1:num_sig1-1
            sig1_all(:,i+1)=sig1(i*split_len+1:1:i*split_len+(split_len),1);
        end

        for i=0:1:num_sig2-1
            sig2_all(:,i+1)=sig2(i*split_len+1:1:i*split_len+(split_len),1);
        end
        
        max_delay = floor(0.8*split_len/10);
        [~,eLag1,eDim1] = phaseSpaceReconstruction(sig1_all,'MaxDim',10,'MaxLag',max_delay);
        [~,eLag2,eDim2] = phaseSpaceReconstruction(sig2_all,'MaxDim',10,'MaxLag',max_delay);


        % traverse through the splited signals
        for k = 1:1:min_sig_size
            sig1_now = sig1_split{k}';
            sig2_now = sig2_split{k}';
            
            % [~,eLag1,eDim1] = phaseSpaceReconstruction(sig1_now');
            % [~,eLag2,eDim2] = phaseSpaceReconstruction(sig2_now');
            sig1_now_embed = delay_embed(sig1_now, eDim1, eLag1);
            sig2_now_embed = delay_embed(sig2_now, eDim2, eLag2);

            sig1_future = delay_embed_future(sig1_now, eLag1);
            sig2_future = delay_embed_future(sig2_now, eLag2);


            TE1 = transfer_entropy(sig1_now_embed,sig2_now_embed,sig1_future,'k',4);
            TE2 = transfer_entropy(sig2_now_embed,sig1_now_embed,sig2_future,'k',4);
            
            % TE1 = transfer_entropy(sig1_now_embed,sig2_now_embed,sig1_future,'xLag', lagSet,'yLag', lagSet, 'wLag', lagSet,'k',4);
            % TE2 = transfer_entropy(sig2_now_embed,sig1_now_embed,sig2_future,'xLag', lagSet,'yLag', lagSet, 'wLag', lagSet,'k',4);
            TE1_all = [TE1_all; TE1];
            TE2_all = [TE2_all; TE2];

        end
        results = [results mean(TE2_all)];
        results = [results mean(TE1_all)];

        ava_ch_index = ava_ch_index+1;
    end
    ava_ch_index = 1;
end
end_index = size(results,2);
idx = [start_index,end_index];