%% Demo code for using the toolbox
clear
clc


%%
signal = rand(10000,5,30);                                                  % assume you have 5 processes (EEG channels), each consist of 1000 time points, and 10 trials (subjects/patients)

epoch_len = 1000;                                                           % here we define 1000 datapoint as 1 epoch
f_s = 250;                                                                  % we assume the sampling frequency is 250 Hz

%% Here we create some filters that will be used in each iteration
delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [12,30];
d_delta = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',delta_bd(1),'CutoffFrequency2',delta_bd(2),'SampleRate',f_s);
d_theta = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',theta_bd(1),'CutoffFrequency2',theta_bd(2),'SampleRate',f_s);
d_alpha = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',alpha_bd(1),'CutoffFrequency2',alpha_bd(2),'SampleRate',f_s);
d_beta = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',beta_bd(1),'CutoffFrequency2',beta_bd(2),'SampleRate',f_s);
filters = {d_delta,d_theta,d_alpha,d_beta};

%%
labels = []; % container for labels
channels = {'Fz';'F3';'C4';'Pz';'P4'};  % define EEG channels

features_AD = [];
features_HC = [];

num_epochs = floor(size(signal,1)/epoch_len);  % we want to further divide the recordings into smaller epochs

for i = 1:size(signal,3)

    % create a label for the subject AD and HC
    if i<=15
        condition = "AD";
        labels = [labels; condition];
    else
        condition = "HC";
        labels = [labels; condition];
    end

    % extract the signal
    signal_subj = signal(:,:,i);

    % ##############   this is THE most important step    #################
    % when formatting the signal, make sure the order is the same as the label
    % Here signal_subj is a 10000 x 5 matrix, we need transform it into a cell

    data_cell = cell(size(signal_subj,2),1);

    for m = 1:num_epochs
        data_cell = cell(size(channels,1),1);
        for n = 1:size(channels,1)
            y = signal_subj(epoch_len*(m-1)+1:epoch_len*m,n)';

            data_cell{n} = double(y)';
        end

        % this is the feature calculation file
        [results, dscrp, idx_ft_all] = Spatial_feature_calculation(data_cell,channels,f_s,filters,[]);

        if condition == "AD"
            features_AD = [features_AD; results];
        elseif condition == "HC"
            features_HC = [features_HC; results];
        end
    end
end

%% Now do classification
sample_num_AD = size(features_AD,1);
sample_num_HC = size(features_HC,1);

min_sample = min([sample_num_HC,sample_num_AD]);

idx_AD = randperm(sample_num_AD,15);                                        % randomly choose a subset of samples
idx_HC = randperm(sample_num_HC,15);

features_AD_rand = features_AD(idx_AD,:);
features_HC_rand = features_HC(idx_HC,:);

feature_space_all = [];
labels_A = string();

for m = 1:size(features_AD_rand,1)
    if m == 1
        labels_A(1) = 'AD';
        labels_A(2) = 'HC';
        feature_space_all = [feature_space_all; features_AD_rand(m,:)];
        feature_space_all = [feature_space_all; features_HC_rand(m,:)];
    else
        labels_A(end+1) = 'AD';
        labels_A(end+1) = 'HC';
        feature_space_all = [feature_space_all; features_AD_rand(m,:)];
        feature_space_all = [feature_space_all; features_HC_rand(m,:)];
    end
end

ft_params = struct();
ft_params.numFeatures = size(dscrp,1);
ft_params.description = dscrp;
ft_params.ID = [1:size(dscrp,1)]';
cfnParams = GiveMeDefaultParams(labels_A);

[ifeat,testStat] = topFeatureHierClust(feature_space_all,cfnParams,ft_params,30); % the last input is find the number of top features

%% Grouped feature classification
group_testStat = groupFeatureHierClust(idx_ft_all,testStat,feature_space_all);

%% Regional feature classification
[dscrp, features_space_all_regional] = Regional_Connectivity(dscrp, feature_space_all);

%%
ft_params = struct();
ft_params.numFeatures = size(dscrp,1);
ft_params.description = dscrp;
ft_params.ID = [1:size(dscrp,1)]';
cfnParams = GiveMeDefaultParams(labels_A);

[ifeat,testStat] = topFeatureHierClust(features_space_all_regional,cfnParams,ft_params,15); % the last input is find the number of top features
