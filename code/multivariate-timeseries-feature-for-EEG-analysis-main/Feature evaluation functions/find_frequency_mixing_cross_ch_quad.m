function [results,dscrp,idx] = find_frequency_mixing_cross_ch_quad(y,ch_labels,split_length,results,dscrp,cfg)
signal = y';
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
size(signal)
signals = {[signal{1} signal{2} signal{3} signal{4} signal{5} signal{6} signal{7} ... 
    signal{8} signal{9} signal{10} signal{11} signal{12} signal{13} signal{14} signal{15} signal{16} signal{17}]};

%% Construct main objects and load data 
% create a data collection and load it with data
dc = freqmix.data.datacollection(signals,...
                                'name','synthetic_scan',... 
                                'channel_names',{'Ch1','Ch2','Ch3','Ch4','Ch5','Ch6','Ch7','Ch8','Ch9','Ch10',...
                                     'Ch11','Ch12','Ch13','Ch14','Ch15','Ch16','Ch17'},...
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


 

quadruplets_HOIs = exp.frequencymixing.quadrupletmixing{1, 1}.hoi;
R1_quadruplets = exp.frequencymixing.quadrupletmixing{1, 1}.R1;
R2_quadruplets = exp.frequencymixing.quadrupletmixing{1, 1}.R2;
E1_quadruplets = exp.frequencymixing.quadrupletmixing{1, 1}.E1;
E2_quadruplets = exp.frequencymixing.quadrupletmixing{1, 1}.E2;


ch1 = exp.frequencymixing.quadrupletmixing{1, 1}.Channel1;
ch2 = exp.frequencymixing.quadrupletmixing{1, 1}.Channel2;
ch3 = exp.frequencymixing.quadrupletmixing{1, 1}.Channel3;
ch4 = exp.frequencymixing.quadrupletmixing{1, 1}.Channel4;

% feature description
for j = 1:size(quadruplets_HOIs,1)
    dscrp{end+1} = string(labels{ch1(j),1})+"_"+string(labels{ch2(j),1})+"_"+string(labels{ch3(j),1})+"_"+...
        string(labels{ch4(j),1})+"_"+num2str(R1_quadruplets(j,1))+"_"+ ... 
        num2str(R2_quadruplets(j,1))+"_"+num2str(E1_quadruplets(j,1))+"_"+num2str(E2_quadruplets(j,1));

    quadruplets_HOI_now = quadruplets_HOIs(j,1);

    results = [results quadruplets_HOI_now];
end



end_index = size(results,2);
idx = [start_index,end_index];