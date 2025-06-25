function [results,dscrp,idx] = find_spectral_TE(y,ch_labels,results,dscrp)
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal);
results = results;
dscrp = dscrp;

TOrder = 6;
FOrder = 6;
TLag = 1;
Type = 1;


start_index = size(results,2)+1;

for m = (1:ava_ch_num-1)
    for n = (m+1:1:ava_ch_num)
        sig1 = signal{m,1};
        sig1_ch = labels{m,1};
        sig2 = signal{n,1};
        sig2_ch = labels{n,1};


        [TransferEntropyM12,~,~,~,~,f] = MainFunctionforFDTES(sig1,sig2,TOrder,FOrder,TLag,Type);

        f = f(1:size(f,1)-FOrder);
        f = round(f);
        % feature description
        TransferEntropyM12 = sum(TransferEntropyM12,2);

        % Find indices where the value of f is between 1 and 4
        delta_index = f >= 1 & f <= 4;

        % Find indices where the value of f is between 5 and 8
        theta_index = f >= 5 & f <= 8;

        % Find indices where the value of f is between 9 and 12
        alpha_index = f >= 9 & f <= 12;

        % Find indices where the value of f is between 13 and 30
        beta_index = f >= 13 & f <= 30;

        TE_delta12 = mean(TransferEntropyM12(delta_index),1);
        TE_theta12 = mean(TransferEntropyM12(theta_index),1);
        TE_alpha12 = mean(TransferEntropyM12(alpha_index),1);
        TE_beta12 = mean(TransferEntropyM12(beta_index),1);


        [TransferEntropyM21,~,~,~,~,~] = MainFunctionforFDTES(sig2,sig1,TOrder,FOrder,TLag,Type);

        f = f(1:size(f,1)-FOrder);
        f = round(f);
        % feature description
        TransferEntropyM21 = sum(TransferEntropyM21,2);

        TE_delta21 = mean(TransferEntropyM21(delta_index),1);
        TE_theta21 = mean(TransferEntropyM21(theta_index),1);
        TE_alpha21 = mean(TransferEntropyM21(alpha_index),1);
        TE_beta21 = mean(TransferEntropyM21(beta_index),1);


        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Spectral_Transfer_Entropy_delta';
        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Spectral_Transfer_Entropy_theta';
        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Spectral_Transfer_Entropy_alpha';
        dscrp{end+1} = string(sig1_ch)+'_to_'+string(sig2_ch)+'_Spectral_Transfer_Entropy_beta';

        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Spectral_Transfer_Entropy_delta';
        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Spectral_Transfer_Entropy_theta';
        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Spectral_Transfer_Entropy_alpha';
        dscrp{end+1} = string(sig2_ch)+'_to_'+string(sig1_ch)+'_Spectral_Transfer_Entropy_beta';

        results = [results TE_delta12 TE_theta12 TE_alpha12 TE_beta12 ...
                    TE_delta21 TE_theta21 TE_alpha21 TE_beta21];

    end
end
end_index = size(results,2);
idx = [start_index,end_index];