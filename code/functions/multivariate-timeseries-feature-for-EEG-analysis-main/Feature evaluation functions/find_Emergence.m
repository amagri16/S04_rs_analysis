function [results,dscrp,idx] = find_Emergence(y, ch_labels, split_length, nplets,filters,results,dscrp,band_name)
 split_len = split_length;
    signal = y;
    labels = ch_labels;
    % how many channels exist
    ava_ch_num = size(signal)-1;

   % we suppose the last input is our marco feature 
    ch_num = 1:ava_ch_num;
    % n-plets to use
    combinations = combnk(ch_num, nplets);
    [rs, ~] = size(combinations);
    start_index = size(results, 2) + 1;

    for i = 1:ava_ch_num-1
        signal{i,1} = filtfilt(filters, signal{i,1});
    end

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

        sigs_ch{end+1} = labels{end, 1};
        
        V_all = y{end};
        V_split = cell(min_sig_size,1);

        % extract the splits for V
        for i = 1:min_sig_size
            V_split{i} = V_all((i-1)*split_len+1:i*split_len, :);
        end

        % Create feature description
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_Psi',band_name];
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_Gamma',band_name];
        dscrp{end+1} = [strjoin(sigs_ch, '_'),'_Delta',band_name];

        % Feature calculation
        Psi_all = [];
        Gamma_all = [];
        Delta_all = [];

        for k = 1:min_sig_size
            sig_now = cellfun(@(x) x{k}, sigs_split, 'UniformOutput', false);
            X = cat(2, sig_now{:});
            V = V_split{k};

           [psi, ~, ~] = EmergencePsi(X, V);
           [gamma] = EmergenceGamma(X, V);
           [delta, ~, ~] = EmergenceDelta(X, V);

            Psi_all = [Psi_all, psi];
            Gamma_all = [Gamma_all, gamma];
            Delta_all = [Delta_all, delta];
        end

        results = [results, mean(Psi_all)];
        results = [results, mean(Gamma_all)];
        results = [results, mean(Delta_all)];
    end
    end_index = size(results, 2);
    idx = [start_index, end_index];
end