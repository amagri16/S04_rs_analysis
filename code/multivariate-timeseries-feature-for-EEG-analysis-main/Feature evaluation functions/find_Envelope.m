function [results, dscrp, idx] = find_Envelope(y,ch_labels,results,dscrp,filters)

signal = y;
% how many channels exist
ava_ch_num = size(signal);
start_index = size(results,2)+1;

for m = 2

    % Band frequency specification
    y = signal{m,1};
    y_33Hz = filtfilt(filters,y);
    y_ch_label = ch_labels{m,1};

    [yupper,ylower] = envelope(y);
    [yupper_33,ylower_33] = envelope(y_33Hz);
    
    pks = findpeaks(y_33Hz);

    results = [results mean(yupper)-mean(ylower)];
    results = [results mean(yupper_33)-mean(ylower_33)];
    results = [results mean(pks)];

    dscrp{end+1} = string(y_ch_label)+'_envelope';
    dscrp{end+1} = string(y_ch_label)+'_33Hz_envelope';
    dscrp{end+1} = string(y_ch_label)+'_33Hz_mean_peak';
end

end_index = size(results,2);
idx = [start_index,end_index];