function [results, dscrp] = corr_coh_test(y,ch_labels,discard,f)
f_s = f;
split_len = 512;
nfft = 0:0.5:30;
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [13,30];
signal = y;
labels = ch_labels;
discard_info = discard;
ava_ch_num = size(signal);
ava_ch_index = 1;
results = [];
dscrp = {};

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};
        sig_size = size(sig1,1);

        % assign the description of the features 
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_pearson_coeff';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_xcorr';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_max_xcorr_value';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_max_xcorr_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_half_max_xorr_lag_width';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_xcov';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_max_xcov_value';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_max_xcov_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_half_max_xcov_lag_width';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_delta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_theta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_alpha_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_mean_beta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_integral_delta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_integral_theta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_integral_alpha_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_global_integral_beta_coh';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_delta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_theta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_alpha_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_beta_coh';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_integral_delta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_integral_theta_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_integral_alpha_coh';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_integral_beta_coh';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_xcorr';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_max_xorr_value';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_max_xorr_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_half_max_xorr_lag_width';

        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_xcov';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_max_xov_value';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_max_xov_lag';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_half_max_xcov_lag_width';
        dscrp{end+1} = string(sig1_ch)+'_vs_'+string(sig2_ch)+'_mean_pearson';

        if discard_info(m) | discard_info(n)
            results = [results nan nan nan nan nan nan nan nan nan nan nan...
                        nan nan nan nan nan nan nan nan nan nan nan...
                        nan nan nan nan nan nan nan nan nan nan nan nan];
        else
            
            % global analysis

            % pearson
            R_pearson_global = corrcoef(sig1,sig2);
            results = [results R_pearson_global(2,1)];

            % global xcorr 
            R = xcorr(sig1,sig2);
            [corr_max,corr_max_index] = max(R);
            max_lag = corr_max_index-sig_size-1;
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);

            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);

            results = [results mean(R) corr_max max_lag half_lag_len];

            % global xcov  
            R = xcov(sig1,sig2);
            [corr_max,corr_max_index] = max(R);
            max_lag = corr_max_index-sig_size-1;
            half_max = 0.5*corr_max;
            % all the index that is larger than half of the max
            corr_half_max_lags = (R>half_max);

            [labeledX, ~] = bwlabel(corr_half_max_lags);
            props = regionprops(labeledX, 'Area', 'PixelList');
            regionLengths = [props.Area];
            half_lag = find(labeledX==labeledX(corr_max_index));
            half_lag_len = size(half_lag,1);

            results = [results mean(R) corr_max max_lag half_lag_len];

            % global coherence



            % global mean bands

            % global integral bands
            cxy = mscohere(sig1,sig2,[],[],nfft,f_s);
            delta_idx = (nfft>delta_bd(1) & nfft<delta_bd(2));
            theta_idx = (nfft>theta_bd(1) & nfft<theta_bd(2));
            alpha_idx = (nfft>alpha_bd(1) & nfft<alpha_bd(2));
            beta_idx = (nfft>beta_bd(1) & nfft<beta_bd(2));

            mean_delta_coh = mean(cxy(delta_idx));
            mean_theta_coh = mean(cxy(theta_idx));
            mean_alpha_coh = mean(cxy(alpha_idx));
            mean_beta_coh = mean(cxy(beta_idx));
                
            integral_delta_coh = sum(cxy(delta_idx));
            integral_theta_coh = sum(cxy(theta_idx));
            integral_alpha_coh = sum(cxy(alpha_idx));
            integral_beta_coh = sum(cxy(beta_idx));

            
            results = [results, mean_delta_coh, mean_theta_coh...
                            mean_alpha_coh, mean_beta_coh...
                            integral_delta_coh, integral_theta_coh...
                            integral_alpha_coh, integral_beta_coh...
                            ];


            sig1_split = {};
            sig2_split = {};
            num_sig1 = floor(size(sig1,1)/split_len);
            num_sig2 = floor(size(sig2,1)/split_len);
            for i=0:1:num_sig1-1
                sig1_split{i+1,1}=sig1(i*split_len+1:1:i*split_len+(split_len));
            end

            for i=0:1:num_sig2-1
                sig2_split{i+1,1}=sig2(i*split_len+1:1:i*split_len+(split_len));
            end

            % xcorr analysis 
            corre_m_n = {};
            corre_m_n_max = {};
            corre_m_n_max_lag = {};
            corre_m_n_half_max_lag_len = {};
            
            % xcov analysis
            cov_m_n = {};
            cov_m_n_max = {};
            cov_m_n_max_lag = {};
            cov_m_n_half_max_lag_len = {};

            % pearson analysis 
            pearson_m_n = {};

            cxy_m_n = {};
           
            corre_m_n_index = 1;
            cxy_m_n_index = 1;
            sig1_split_size = size(sig1_split,1);
            sig2_split_size = size(sig2_split,1);
            min_sig_size = min(sig1_split_size,sig2_split_size);
            % traverse through the splited signals
            
            for k = 1:1:min_sig_size
                % xcor here
                R = xcorr(sig1_split{k,1},sig2_split{k,1});
                [corr_max,corr_max_index] = max(R);
                max_lag = corr_max_index-split_len-1;
                half_max = 0.5*corr_max;
                % all the index that is larger than half of the max
                corr_half_max_lags = (R>half_max);

                [labeledX, ~] = bwlabel(corr_half_max_lags);
                props = regionprops(labeledX, 'Area', 'PixelList');
                regionLengths = [props.Area];
                half_lag = find(labeledX==labeledX(corr_max_index));
                half_lag_len = size(half_lag,1);

                corre_m_n_max{end+1} = corr_max;
                corre_m_n_max_lag{end+1} = max_lag;
                corre_m_n_half_max_lag_len{end+1} = half_lag_len;
            
                corre_m_n{k,1} = mean(R);

                % xcov here
                R = xcov(sig1_split{k,1},sig2_split{k,1});
                [corr_max,corr_max_index] = max(R);
                max_lag = corr_max_index-split_len-1;
                half_max = 0.5*corr_max;
                % all the index that is larger than half of the max
                corr_half_max_lags = (R>half_max);

                [labeledX, ~] = bwlabel(corr_half_max_lags);
                props = regionprops(labeledX, 'Area', 'PixelList');
                regionLengths = [props.Area];
                half_lag = find(labeledX==labeledX(corr_max_index));
                half_lag_len = size(half_lag,1);

                cov_m_n_max{end+1} = corr_max;
                cov_m_n_max_lag{end+1} = max_lag;
                cov_m_n_half_max_lag_len{end+1} = half_lag_len;
            
                cov_m_n{k,1} = mean(R);

                
                R = corrcoef(sig1_split{k,1},sig2_split{k,1});
                pearson_m_n{k,1} = R(2,1);

                corre_m_n_index = corre_m_n_index+1;
                cxy_m_n{k,1} = mscohere(sig1_split{k,1},sig2_split{k,1},[],[],nfft,f_s);
                cxy_m_n_index = cxy_m_n_index+1;
            end
            corr_results = mean(transpose(cell2mat(corre_m_n)));
            cov_results = mean(transpose(cell2mat(cov_m_n)));
            pearson_results = mean(transpose(cell2mat(pearson_m_n)));
            cxy_table = mean(transpose(cell2mat(cxy_m_n)),2);
            %% coherence on each band
            delta_idx = (nfft>delta_bd(1) & nfft<delta_bd(2));
            theta_idx = (nfft>theta_bd(1) & nfft<theta_bd(2));
            alpha_idx = (nfft>alpha_bd(1) & nfft<alpha_bd(2));
            beta_idx = (nfft>beta_bd(1) & nfft<beta_bd(2));

            mean_delta_coh = mean(cxy_table(delta_idx));
            mean_theta_coh = mean(cxy_table(theta_idx));
            mean_alpha_coh = mean(cxy_table(alpha_idx));
            mean_beta_coh = mean(cxy_table(beta_idx));

            integral_delta_coh = sum(cxy_table(delta_idx));
            integral_theta_coh = sum(cxy_table(theta_idx));
            integral_alpha_coh = sum(cxy_table(alpha_idx));
            integral_beta_coh = sum(cxy_table(beta_idx));
                

                results = [results, mean_delta_coh, mean_theta_coh,...
                            mean_alpha_coh, mean_beta_coh ...
                            integral_delta_coh, integral_theta_coh...
                            integral_alpha_coh, integral_beta_coh...
                            corr_results...
                            mean(cell2mat(corre_m_n_max),2)...
                            mean(cell2mat(corre_m_n_max_lag),2)...
                            mean(cell2mat(corre_m_n_half_max_lag_len),2)...
                            cov_results...
                            mean(cell2mat(cov_m_n_max),2)...
                            mean(cell2mat(cov_m_n_max_lag),2)...
                            mean(cell2mat(cov_m_n_half_max_lag_len),2)...
                            pearson_results];
            ava_ch_index = ava_ch_index+1;
        end
        ava_ch_index = 1;
    end
end
dscrp = transpose(dscrp);