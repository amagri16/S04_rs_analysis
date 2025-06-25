function [results, dscrp, idx] = find_O_info(y, ch_labels, split_length, nplets,filters,results,dscrp,band_name)
    split_len = split_length;
    signal = y;
    labels = ch_labels;
    % how many channels exist
    ava_ch_num = size(signal); 
    ch_num = 1:ava_ch_num;
    % n-plets to use
    combinations = combnk(ch_num, nplets);
    [rs, ~] = size(combinations);
    start_index = size(results, 2) + 1;

    % for i = 1:ava_ch_num-1
    %     signal{i,1} = filtfilt(filters, signal{i,1});
    % end

    for m = 1:rs
        sigs = cell(nplets, 1);
        sigs_ch = cell(nplets, 1);
        sigs_split = cell(nplets, 1);
        min_sig_size = inf;

        % Extract signals and labels for each combination
        for n = 1:nplets
            sigs{n} = signal{combinations(m, n), 1};
            sigs_ch{n} = labels{combinations(m, n), 1};

            % Split signals
            num_sigs = floor(size(sigs{n}, 1) / split_len);
            min_sig_size = min(min_sig_size, num_sigs);
            sigs_split{n} = mat2cell(sigs{n}(1:num_sigs*split_len, :), repmat(split_len, num_sigs, 1), size(sigs{n}, 2));
        end

        % Create feature description
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_O_Information',band_name];
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_S_Information',band_name];
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_Total_correlation',band_name];
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_Dual_total_correlation',band_name];

        % Feature calculation
        OI_all = [];
        SI_all = [];
        I_all = [];
        D_all = [];

        for k = 1:min_sig_size
            sig_now = cellfun(@(x) x{k}, sigs_split, 'UniformOutput', false);
            [O, I, D] = calcO_logdet(cov(cat(2, sig_now{:})));
            OI_all = [OI_all, O];
            SI_all = [SI_all, I+D];
            I_all = [I_all, I];
            D_all = [D_all D];
        end

        results = [results, mean(OI_all)];
        results = [results, mean(SI_all)];
        results = [results, mean(I_all)];
        results = [results, mean(D_all)];
    end
    end_index = size(results, 2);
    idx = [start_index, end_index];
end
