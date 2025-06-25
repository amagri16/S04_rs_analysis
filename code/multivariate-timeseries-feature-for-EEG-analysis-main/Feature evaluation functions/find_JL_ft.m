function [results, dscrp, idx] = find_JL_ft(y, Fs, ch_labels, results, dscrp, filters)
    % find_JL_ft: extract global EEG features per channel and build descriptions
    fprintf('*** Running find_JL_ft from: %s ***\n', mfilename('fullpath'));
    start_index = size(results,2) + 1;
    nCh = numel(y);

    for ch = 1:nCh
        sig = y{ch};
        if size(sig,2) < size(sig,1)
            sig = sig';
        end
        lbl = ch_labels{ch};

        %% --- PSD & band powers ---
        [psd_y, f] = pmtm(sig, 4, length(sig), Fs);
        psd_norm = psd_y / sum(psd_y);
        delta_m = f>1   & f<4;
        theta_m = f>4   & f<8;
        alpha_m = f>8   & f<12;
        beta_m  = f>13  & f<30;
        b_n_m   = f>25  & f<29;
        gamma_m = f>30  & f<35;
        dp   = mean(psd_norm(delta_m), 'omitnan');
        tp   = mean(psd_norm(theta_m), 'omitnan');
        ap   = mean(psd_norm(alpha_m), 'omitnan');
        bp   = mean(psd_norm(beta_m),  'omitnan');
        bnp  = mean(psd_norm(b_n_m),   'omitnan');
        gp   = mean(psd_norm(gamma_m), 'omitnan');
        results = [results, log([dp,tp,ap,bp,bnp,gp]), ...
                   log([dp/tp, dp/ap, dp/bp, tp/ap, tp/bp, ap/bp, (dp+tp)/(ap+bp)])];
        masks = {delta_m, theta_m, alpha_m, beta_m, b_n_m, gamma_m};
        for msk = masks
            v = psd_norm(msk{1}); fr = f(msk{1}); [~,iMax] = max(v);
            results(end+1) = fr(iMax);
        end

        %% --- Temporal coherence ---
        epSec = 1; nEpoch = floor(length(sig)/Fs/epSec) - 1; C = []; count = 0;
        for i=1:nEpoch-1
            for j=i+1:nEpoch
                count = count + 1;
                seg1 = sig((i-1)*Fs*epSec+1:i*Fs*epSec);
                seg2 = sig((j-1)*Fs*epSec+1:j*Fs*epSec);
                [C(count,:), fC] = mscohere(seg1, seg2, [], [], 0:0.5:30, Fs);
            end
        end
        Cmean = mean(C, 'omitnan');
        results = [results, mean(Cmean(fC>1  & fC<4 ), 'omitnan'), ...
                          mean(Cmean(fC>4  & fC<8 ), 'omitnan'), ...
                          mean(Cmean(fC>8  & fC<12),'omitnan'), ...
                          mean(Cmean(fC>13 & fC<30),'omitnan')];

        %% --- Phase-Amplitude Coupling (PAC) ---
        d_sig = filtfilt(filters{1}, 1, sig);
        t_sig = filtfilt(filters{2}, 1, sig);
        a_sig = filtfilt(filters{3}, 1, sig);
        b_sig = filtfilt(filters{4}, 1, sig);
        ph_d  = angle(hilbert(d_sig));  amp_t = abs(hilbert(t_sig));
        ph_t  = angle(hilbert(t_sig));  amp_a = abs(hilbert(a_sig));
        ph_a  = angle(hilbert(a_sig));  amp_b = abs(hilbert(b_sig));
        pairs = {{ph_d,amp_t},{ph_d,amp_a},{ph_d,amp_b},{ph_t,amp_a},{ph_t,amp_b},{ph_a,amp_b}};
        pac = zeros(1, numel(pairs));
        for p=1:numel(pairs)
            ph = pairs{p}{1}; am = pairs{p}{2}; idxC = find_cycle(ph);
            zt = am(1:idxC) .* exp(1j*ph(1:idxC));
            pac(p) = abs(sum(zt)/sum(am(1:idxC)));
        end
        results = [results, pac];

        %% --- (Optional) Catch-22 features ---
        try
            sig = double(sig);
            [cV, cN] = catch22_all(sig);
            results = [results, cV'];
            pfx = lbl + "_";
            dscrp = [dscrp, strcat(pfx, cN')'];
        catch ME
            warning('%s: catch22_all unavailable; skipping those features. %s', ME.identifier, ME.message);
        end

        %% --- Build feature names (per channel) ---
        pfx = lbl + "_";
        dscrp = [ dscrp, ...
            pfx+"delta", pfx+"theta", pfx+"alpha", pfx+"beta", pfx+"beta_narrow", pfx+"gamma", ...
            pfx+"DR", pfx+"AR", pfx+"BR", pfx+"TR", pfx+"TBR", pfx+"ABR", pfx+"TransR", ...
            pfx+"peak_delta", pfx+"peak_theta", pfx+"peak_alpha", pfx+"peak_beta", pfx+"peak_beta_narrow", pfx+"peak_gamma", ...
            pfx+"TCdelta", pfx+"TCtheta", pfx+"TCalpha", pfx+"TCbeta", ...
            pfx+"PAC_DT", pfx+"PAC_DA", pfx+"PAC_DB", pfx+"PAC_TA", pfx+"PAC_TB", pfx+"PAC_AB" ...
        ];
    end
    idx = [start_index, size(results,2)];
end

%% Helper: find first zero-crossing cycle
function idx = find_cycle(ph)
    base = ph(1);
    for k = 2:numel(ph)
        if (ph(k-1)>base && ph(k)<=base) || (ph(k-1)<base && ph(k)>=base)
            idx = k; return;
        end
    end
    idx = numel(ph);
end