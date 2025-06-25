function [results,dscrp,idx] = find_frequency_mixing_single_ch(y,ch_labels,split_length,results,dscrp,cfg)
signal = y;
labels = ch_labels;
% how many channels exist
ava_ch_num = size(signal,1);
results = results;
dscrp = dscrp;

start_index = size(results,2)+1;


% set hyperparameters for data
fs = 256; % sampling frequency
time = round(split_length/fs);
n_channels = 1; % number of channels
n_samples = size(signal); % number of channels here. This might be a bit confusing but each sample is one channel, easier to extract info


patient_ids = 1:ava_ch_num;
signal_ids = 1:ava_ch_num;
group_labels = ones(1,ava_ch_num);
signal_names = labels';
signals = signal';

%% Construct main objects and load data 
% create a data collection and load it with data
dc = freqmix.data.datacollection(signals,...
                                'name','synthetic_scan',... 
                                'channel_names',{'Ch'},...
                                'group_id', group_labels,...
                                'signal_ids', signal_ids,...
                                'sample_ids', patient_ids);

% construct experiment with configuration
exp = freqmix.experiments.experiment('name','synthetic_data',...
                                     'config',cfg,...
                                     'datacollection',dc);

%% run main frequency mixing functions                                 
% initiate frequency mixing computations

exp = exp.run_frequencymixing();

%%
for i = 1:ava_ch_num
    sig1_ch = labels{i,1};

    triplets_HOIs = exp.frequencymixing.tripletmixing{i, 1}.hoi;
    harmonics_HOIs = exp.frequencymixing.harmonicmixing{i, 1}.hoi;

    Freq1_triplets = exp.frequencymixing.tripletmixing{i, 1}.Freq1;
    Freq2_triplets = exp.frequencymixing.tripletmixing{i, 1}.Freq2;
    Freq3_triplets = exp.frequencymixing.tripletmixing{i, 1}.Freq3;

    Freq1_harmonics = exp.frequencymixing.harmonicmixing{i, 1}.Freq1;
    Freq2_harmonics = exp.frequencymixing.harmonicmixing{i, 1}.Freq2;

    % feature description
    for j = 1:size(triplets_HOIs,1)
        dscrp{end+1} = string(sig1_ch)+"_frequency_mixing_triplets_"+...
            num2str(Freq1_triplets(j,1))+"_"+num2str(Freq2_triplets(j,1))+"_"+num2str(Freq3_triplets(j,1));
        
        triplets_HOI_now = triplets_HOIs(j,1);

        results = [results triplets_HOI_now];
    end

    for k = 1:size(harmonics_HOIs,1)
        dscrp{end+1} = string(sig1_ch)+"_frequency_mixing_harmonics_"+...
            num2str(Freq1_harmonics(k,1))+"_"+num2str(Freq2_harmonics(k,1));

        Harmonics_HOIs_now = harmonics_HOIs(k,1);

        results = [results Harmonics_HOIs_now];
    end

end


end_index = size(results,2);
idx = [start_index,end_index];