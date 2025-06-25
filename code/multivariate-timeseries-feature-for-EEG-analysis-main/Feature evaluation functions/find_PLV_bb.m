function [results,dscrp,idx] = find_PLV_bb(y,ch_labels,split_length,filters,results,dscrp)
% Author: Tianyu Wei/Junheng Li
split_len = split_length;
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
results = results;
dscrp = dscrp;

f_interest = [2,4,8,12,16,20,24,28,30];
combinations = nchoosek(f_interest, 2);

start_index = size(results,2)+1;

for k = 1:size(combinations,1)
    for m = 1:ava_ch_num
        for n = 1:ava_ch_num
            if m~=n
                s1 = signal{m,1};
                s1 = s1(2:end) - s1(1:end-1);
                s1 = [s1;0];

                s2 = signal{n,1};
                s2 = s2(2:end) - s2(1:end-1);
                s2 = [s2;0];

                sig1 = filtfilt(filters{combinations(k,1)}, s1);
                sig1_ch = labels{m,1};
                sig2 = filtfilt(filters{combinations(k,2)}, s2);
                sig2_ch = labels{n,1};

                difff_upper = combinations(k,2)-combinations(k,1)+4+1;
                difff_lower = combinations(k,2)-combinations(k,1)-4;

                if difff_lower <= 0
                    difff_lower = 0;
                end
                
                difff_lower = difff_lower+1;

                summ_upper = combinations(k,2)+combinations(k,1)+4+1;
                summ_lower = combinations(k,2)+combinations(k,1)-4+1;

                sig1_diff = filtfilt(filters{1,difff_upper}, signal{m,1});
                sig1_diff = filtfilt(filters{2,difff_lower}, sig1_diff);
                sig1_sum = filtfilt(filters{1,summ_upper}, signal{m,1});
                sig1_sum = filtfilt(filters{2,summ_lower}, sig1_sum);

                sig2_diff = filtfilt(filters{1,difff_upper}, signal{n,1});
                sig2_diff = filtfilt(filters{2,difff_lower}, sig1_diff);
                sig2_sum = filtfilt(filters{1,summ_upper}, signal{n,1});
                sig2_sum = filtfilt(filters{2,summ_lower}, sig1_sum);


                % feature description
                dscrp{end+1} = string(sig1_ch)+'_'+num2str(combinations(k,1))+ ...
                    '_Hz'+'_vs_'+string(sig2_ch)+'_'+num2str(combinations(k,2))+'_Hz'+'_PLV_'+'diff_'+string(sig1_ch);

                dscrp{end+1} = string(sig1_ch)+'_'+num2str(combinations(k,1))+ ...
                    '_Hz'+'_vs_'+string(sig2_ch)+'_'+num2str(combinations(k,2))+'_Hz'+'_PLV_'+'sum_'+string(sig1_ch);

                dscrp{end+1} = string(sig1_ch)+'_'+num2str(combinations(k,1))+ ...
                    '_Hz'+'_vs_'+string(sig2_ch)+'_'+num2str(combinations(k,2))+'_Hz'+'_PLV_'+'diff_'+string(sig2_ch);

                dscrp{end+1} = string(sig1_ch)+'_'+num2str(combinations(k,1))+ ...
                    '_Hz'+'_vs_'+string(sig2_ch)+'_'+num2str(combinations(k,2))+'_Hz'+'_PLV_'+'sum_'+string(sig2_ch);
                
                
                % feature calculation container for each paired-channels
                PLV_all_diff_sig1 = [];
                PLV_all_sum_sig1 = [];
                PLV_all_diff_sig2 = [];
                PLV_all_sum_sig2 = [];

                sig1_split = {};
                sig2_split = {};

                sig1_diff_split = {};
                sig1_sum_split = {};
                sig2_diff_split = {};
                sig2_sum_split = {};

                num_sig1 = floor(size(sig1,1)/split_len);
                num_sig2 = floor(size(sig2,1)/split_len);

                % split the data into sub parts
                for i=0:1:num_sig1-1
                    sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
                    sig1_diff_split{i+1,1}=sig1_diff(i*split_len+1:1:i*split_len+(split_len));
                    sig1_sum_split{i+1,1}=sig1_sum(i*split_len+1:1:i*split_len+(split_len));
                end

                for i=0:1:num_sig2-1
                    sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
                    sig2_diff_split{i+1,1}=sig2_diff(i*split_len+1:1:i*split_len+(split_len));
                    sig2_sum_split{i+1,1}=sig2_sum(i*split_len+1:1:i*split_len+(split_len));
                end

                sig1_split_size = size(sig1_split,1);
                sig2_split_size = size(sig2_split,1);
                min_sig_size = min(sig1_split_size,sig2_split_size);
                % traverse through the splited signals

                for j = 1:1:min_sig_size
                    % add feature analysis for each signal pearson
                    sig1_now = sig1_split{j}';
                    sig2_now = sig2_split{j}';

                    sig1_diff_now = sig1_diff_split{j}';
                    sig1_sum_now = sig1_sum_split{j}';
                    sig2_diff_now = sig2_diff_split{j}';
                    sig2_sum_now = sig2_sum_split{j}';

                    % block1 
                    phase_sig1 = angle(hilbert(sig1_now));
                    phase_sig2 = angle(hilbert(sig2_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_a = abs(sum(e,2)) / split_len;
                   
                    phase_sig1 = angle(hilbert(sig1_now));
                    phase_sig2 = angle(hilbert(sig1_diff_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_b = abs(sum(e,2)) / split_len;

                    phase_sig1 = angle(hilbert(sig2_now));
                    phase_sig2 = angle(hilbert(sig1_diff_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_c = abs(sum(e,2)) / split_len;

                    PLV_diff_sig1 = mean([plv_a,plv_b,plv_c]);

                    % block2
                    phase_sig1 = angle(hilbert(sig1_now));
                    phase_sig2 = angle(hilbert(sig1_sum_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_b = abs(sum(e,2)) / split_len;

                    phase_sig1 = angle(hilbert(sig2_now));
                    phase_sig2 = angle(hilbert(sig1_sum_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_c = abs(sum(e,2)) / split_len;

                    PLV_sum_sig1 = mean([plv_a,plv_b,plv_c]);

                    % block3
                    phase_sig1 = angle(hilbert(sig1_now));
                    phase_sig2 = angle(hilbert(sig2_diff_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_b = abs(sum(e,2)) / split_len;

                    phase_sig1 = angle(hilbert(sig2_now));
                    phase_sig2 = angle(hilbert(sig2_diff_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_c = abs(sum(e,2)) / split_len;

                    PLV_diff_sig2 = mean([plv_a,plv_b,plv_c]);

                    % block 4
                    phase_sig1 = angle(hilbert(sig1_now));
                    phase_sig2 = angle(hilbert(sig2_sum_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_b = abs(sum(e,2)) / split_len;

                    phase_sig1 = angle(hilbert(sig2_now));
                    phase_sig2 = angle(hilbert(sig2_sum_now));
                    e = exp(1i*(phase_sig1 - phase_sig2));
                    plv_c = abs(sum(e,2)) / split_len;

                    PLV_sum_sig2 = mean([plv_a,plv_b,plv_c]);


                    % block finished 
                    PLV_all_diff_sig1 = [PLV_all_diff_sig1 PLV_diff_sig1];
                    PLV_all_sum_sig1 = [PLV_all_diff_sig1 PLV_sum_sig1];
                    PLV_all_diff_sig2 = [PLV_all_diff_sig1 PLV_diff_sig2];
                    PLV_all_sum_sig2 = [PLV_all_diff_sig1 PLV_sum_sig2];
                end
                results = [results mean(PLV_all_diff_sig1) mean(PLV_all_sum_sig1) ...
                    mean(PLV_all_diff_sig2) mean(PLV_all_sum_sig2)];
            end
        end
    end
end
end_index = size(results,2);
idx = [start_index,end_index];