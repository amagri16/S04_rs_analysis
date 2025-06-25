%% master_eeg_analysis_pipeline.m
% Complete automated pipeline: feature extraction → selection → model evaluation
% for all cognitive tasks, with final bar plot and feature selection table

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions', 'multivariate-timeseries-feature-for-EEG-analysis-main', 'HCTSA classification'));
disp(which('topFeatureHierClust'));

clear; clc;

%% 0) Setup and initialization
fprintf('=== Master EEG Analysis Pipeline ===\n');

% Base directories
baseDataDir = '/Users/andreamagri/Desktop/thesis/S04_analysis/data/processed';
[scriptPath,~,~] = fileparts(mfilename('fullpath'));
rootDir = fullfile(scriptPath,'..','..');

% Add required paths
fEvalDir = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','Feature evaluation functions');
fHCTSADir = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','HCTSA classification');
if exist(fEvalDir,'dir'), addpath(genpath(fEvalDir)); end
if exist(fHCTSADir,'dir'), addpath(genpath(fHCTSADir)); end
addpath(fullfile(baseDataDir,'..','..','code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','Feature evaluation functions'));
addpath(fullfile(rootDir, 'code', 'functions', 'multivariate-timeseries-feature-for-EEG-analysis-main', 'HCTSA classification'));
addpath(genpath('/Users/andreamagri/Desktop/thesis/S04_analysis/code/functions/multivariate-timeseries-feature-for-EEG-analysis-main/Feature evaluation functions/catch22-master/wrap_Matlab'))

% Create output directories
featDir = fullfile(rootDir,'code','data','features');
if ~exist(featDir,'dir'), mkdir(featDir); end

%% 1) Define all performance series (from performance_feature_extraction.m)
% Raw z-scored 10-session vectors
om   = [ 0.3825, -0.7015,  0.3825, -0.0752,  1.7555,  -0.8889,  0.8402,  0.7438, -0.9906, -1.4483];
omd  = [ 0.3753, -1.0148,  0.9546, -0.8063,  0.0741,  -0.5514,  0.9546,  1.8350, -1.0148, -0.8063];
wm   = [ 0.3446, -1.4165,  0.7848, -0.2106, -1.5888,   0.7083, -0.0957,  0.3446, -0.5360,  1.6654];
wmd  = [ 0.5697, -2.1105,  0.0208,  0.0208,  0.7312,   0.7312, -0.3344,  1.3953, -1.0449,  0.0208];
fnat = [-1.5479, -0.6669,  1.1391, -0.6669, -1.2689,   1.2496,  0.5371,  0.9548,  0.0702,  0.1999];
fnatd= [ 0.8568,  0.0261,  0.1100,  1.1410, -1.3365,   0.8712,  0.2734,  0.1100,  0.0261, -2.0781];
mst  = [ 1.5571, -0.4104,  0.7878,  0.0546,  0.5566,  -1.3808, -0.6293, -0.5994, -1.1534, -1.4483];
pca1 = [-0.6360,  3.0049, -1.6104,  0.7023, -0.0384,  -0.3482, -1.0221, -2.8583,  2.0115,  0.7947];

% Pack all series
scores = struct( ...
    'pca1',              pca1, ...
    'om',                om, ...
    'omd',               omd, ...
    'wm',                wm, ...
    'wmd',               wmd, ...
    'fnat',              fnat, ...
    'fnatd',             fnatd, ...
    'mst',               mst);

% Compute composites
scores.composite_om       = (scores.om   + scores.omd)   / 2;
scores.composite_wm       = (scores.wm   + scores.wmd)   / 2;
scores.composite_fnat     = (scores.fnat + scores.fnatd) / 2;
scores.composite_immediate = (scores.om + scores.wm + scores.fnat) / 3;
scores.composite_delayed   = (scores.omd + scores.wmd + scores.fnatd) / 3;

% Get all task names
allTasks = fieldnames(scores);
nTasks = numel(allTasks);

fprintf('Processing %d tasks: %s\n', nTasks, strjoin(allTasks, ', '));

%% 2) Initialize parallel pool for feature extraction
try
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    use_parallel = true;
catch
    warning('Parallel toolbox unavailable; using sequential mode.');
    use_parallel = false;
end

%% 3) Feature extraction parameters (from performance_feature_extraction.m)
Fs_target = 1000;
epoch_length = 5 * Fs_target;
bands = {[1 4],[4 8],[8 12],[13 30]};
order = 64;
filters = cell(1,numel(bands));
for i = 1:numel(bands)
    f_nyq = Fs_target/2;
    bn = bands{i}/f_nyq;
    filters{i} = fir1(order, bn, 'bandpass', hamming(order+1));
end

%% 4) Feature selection parameters (from regional_uncorrelated_top10.m)
corrThresh = 0.50;    % maximum allowed abs(corr) between any two kept features
accThresh = 61;      % only keep features with accuracy ≥ this percent
maxKeep = 10;        % stop after this many survivors

%% 5) Model evaluation parameters (from train_and_evaluate_models.m)
kFold = 5;
rng(1,'twister');
models = {'LDA','Logistic','SVM'};

%% 6) Initialize result storage
allAccuracies = struct();
allSelectedFeatures = struct();
session_types = {'pre_day','stim_day3','post_day','follow_up1','follow_up3'};

%% 7) Main processing loop
fprintf('\n=== Processing Tasks ===\n');
for taskIdx = 1:nTasks
    selName = allTasks{taskIdx};
    fprintf('\n--- Processing task: %s (%d/%d) ---\n', selName, taskIdx, nTasks);
    
    %% STEP 1: Feature Extraction (adapted from performance_feature_extraction.m)
    fprintf('Step 1: Feature extraction...\n');
    
    % Get performance scores for this task
    all_scores = scores.(selName);
    
    % Create ±0.5-z labels
    z_hi = 0.25;
    z_lo = -0.25;
    is_high = all_scores >= z_hi;
    is_low = all_scores <= z_lo;
    keep = is_high | is_low;
    
    perf_lbl = is_high(keep);
    all_scores = all_scores(keep);
    
    % Build session_info
    session_info = {};
    idx = 1;
    lab_ptr = 1;
    
    for w = 1:2
        for d = 1:numel(session_types)
            if ~keep(idx)
                idx = idx + 1;
                continue
            end
            session_info{end+1,1} = struct( ...
                'week', sprintf('week_%d',w), ...
                'day', session_types{d}, ...
                'score', all_scores(lab_ptr), ...
                'is_high_performance', perf_lbl(lab_ptr) );
            idx = idx + 1;
            lab_ptr = lab_ptr + 1;
        end
    end
    
    % Process all sessions for this task
    N = numel(session_info);
    session_feats = cell(N,1);
    session_labels = cell(N,1);
    session_meta = cell(N,1);
    name_list = cell(N,1);
    
    if use_parallel
        parfor i = 1:N
            [F,L,mi,nl] = process_session(session_info{i}, Fs_target, epoch_length, filters, baseDataDir);
            session_feats{i} = F;
            session_labels{i} = L;
            session_meta{i} = mi;
            name_list{i} = nl;
        end
    else
        for i = 1:N
            [F,L,mi,nl] = process_session(session_info{i}, Fs_target, epoch_length, filters, baseDataDir);
            session_feats{i} = F;
            session_labels{i} = L;
            session_meta{i} = mi;
            name_list{i} = nl;
        end
    end
    
    % Concatenate results
    gthr = ~cellfun(@isempty, session_feats);
    if ~any(gthr)
        warning('No valid sessions found for task %s. Skipping.', selName);
        continue;
    end
    
    allF = vertcat(session_feats{gthr});
    allLbls = vertcat(session_labels{gthr});
    allMeta = vertcat(session_meta{gthr});
    ds = name_list{find(~cellfun(@isempty,name_list),1)};
    
    high_feats = allF(allLbls==1, :);
    low_feats = allF(allLbls==0, :);
    
    % Check if we have enough data
    if size(high_feats,1) < 2 || size(low_feats,1) < 2
        warning('Insufficient data for task %s (high: %d, low: %d). Skipping.', ...
            selName, size(high_feats,1), size(low_feats,1));
        continue;
    end
    
    fprintf('  High: %d epochs, Low: %d epochs, Features: %d\n', ...
        size(high_feats,1), size(low_feats,1), numel(ds));
    
    %% STEP 2: Feature Selection (adapted from regional_uncorrelated_top10.m)
    fprintf('Step 2: Feature selection...\n');
    
    % Assemble data for feature ranking
    X = [high_feats; low_feats];
    Y = [ones(size(high_feats,1),1); zeros(size(low_feats,1),1)];
    numFeat = size(X,2);
    
    % SVM ranking parameters
    cfnParams.numClasses = 2;
    cfnParams.cat_labels = Y;
    cfnParams.whatLossFun = 'classificationError';
    cfnParams.whatLossUnits = '%';
    
    fparams.numFeatures = numFeat;
    fparams.ID = (1:numFeat)';
    fparams.description = ds;
    
    % Rank all features via single-feature linear SVM
    topN = 10;
    [topIdx, perFeatureAcc] = topFeatureHierClust(X, cfnParams, fparams, topN);
    
    % Greedy uncorrelated + high-accuracy selection
    selectedIdx = [];
    for ii = 1:numFeat
        f = topIdx(ii);
        acc = perFeatureAcc(f);
        if acc < accThresh
            break;
        end
        if isempty(selectedIdx)
            selectedIdx(end+1) = f;
        else
            r = abs(corr(X(:,f), X(:,selectedIdx)));
            if all(r < corrThresh)
                selectedIdx(end+1) = f;
            end
        end
        if numel(selectedIdx) == maxKeep
            break;
        end
    end
    
    % Ensure minimum 5 features
    if numel(selectedIdx) < 5
        remaining = setdiff(topIdx, selectedIdx, 'stable');
        nToAdd = 5 - numel(selectedIdx);
        toAdd = remaining(1:min(nToAdd, numel(remaining)));
        selectedIdx = [selectedIdx(:); toAdd(:)];
    end
    
    finalIdx = selectedIdx;
    finalNames = ds(finalIdx);
    finalAcc = perFeatureAcc(finalIdx);
    
    fprintf('  Selected %d features\n', numel(finalIdx));
    
    %% STEP 3: Model Training & Evaluation (adapted from train_and_evaluate_models.m)
    fprintf('Step 3: Model evaluation...\n');
    
    X_sel = X(:, finalIdx);
    
    % Cross-validation setup
    cvp = cvpartition(Y,'KFold',kFold);
    
    % Initialize metrics storage
    acc = struct();
    AUC = struct();
    precision = struct();
    recall = struct();
    f1score = struct();
    
    % Cross-validation loop
    for fold = 1:kFold
        tr = training(cvp,fold);
        te = test(cvp,fold);
        Xtr = X_sel(tr,:); Ytr = Y(tr);
        Xte = X_sel(te,:); Yte = Y(te);
        
        % LDA
        try
            mdlLDA = fitcdiscr(Xtr, Ytr);
            [predLDA, scoreLDA] = predict(mdlLDA, Xte);
            acc.LDA(fold) = mean(predLDA==Yte);
            confMat = confusionmat(Yte,predLDA);
            [~,~,~,AUC.LDA(fold)] = perfcurve(Yte,scoreLDA(:,2),1);
            if size(confMat,1) == 2 && size(confMat,2) == 2
                TP = confMat(2,2); FP = confMat(1,2); FN = confMat(2,1);
                precision.LDA(fold) = TP/(TP+FP);
                recall.LDA(fold) = TP/(TP+FN);
                f1score.LDA(fold) = 2*(precision.LDA(fold)*recall.LDA(fold))/(precision.LDA(fold)+recall.LDA(fold));
            end
        catch
            acc.LDA(fold) = NaN;
            AUC.LDA(fold) = NaN;
            precision.LDA(fold) = NaN;
            recall.LDA(fold) = NaN;
            f1score.LDA(fold) = NaN;
        end
        
        % Logistic Regression
        try
            mdlLog = fitglm(Xtr, Ytr, 'Distribution','binomial','Link','logit');
            scoresLog = predict(mdlLog, Xte);
            predLog = double(scoresLog>=0.5);
            acc.Logistic(fold) = mean(predLog==Yte);
            confMat = confusionmat(Yte,predLog);
            [~,~,~,AUC.Logistic(fold)] = perfcurve(Yte,scoresLog,1);
            if size(confMat,1) == 2 && size(confMat,2) == 2
                TP = confMat(2,2); FP = confMat(1,2); FN = confMat(2,1);
                precision.Logistic(fold) = TP/(TP+FP);
                recall.Logistic(fold) = TP/(TP+FN);
                f1score.Logistic(fold) = 2*(precision.Logistic(fold)*recall.Logistic(fold))/(precision.Logistic(fold)+recall.Logistic(fold));
            end
        catch
            acc.Logistic(fold) = NaN;
            AUC.Logistic(fold) = NaN;
            precision.Logistic(fold) = NaN;
            recall.Logistic(fold) = NaN;
            f1score.Logistic(fold) = NaN;
        end
        
        % SVM
        try
            mdlSVM = fitcsvm(Xtr, Ytr, 'KernelFunction','linear','Standardize',true);
            [predSVM, scoreSVM] = predict(mdlSVM, Xte);
            acc.SVM(fold) = mean(predSVM==Yte);
            confMat = confusionmat(Yte,predSVM);
            [~,~,~,AUC.SVM(fold)] = perfcurve(Yte,scoreSVM(:,2),1);
            if size(confMat,1) == 2 && size(confMat,2) == 2
                TP = confMat(2,2); FP = confMat(1,2); FN = confMat(2,1);
                precision.SVM(fold) = TP/(TP+FP);
                recall.SVM(fold) = TP/(TP+FN);
                f1score.SVM(fold) = 2*(precision.SVM(fold)*recall.SVM(fold))/(precision.SVM(fold)+recall.SVM(fold));
            end
        catch
            acc.SVM(fold) = NaN;
            AUC.SVM(fold) = NaN;
            precision.SVM(fold) = NaN;
            recall.SVM(fold) = NaN;
            f1score.SVM(fold) = NaN;
        end
    end
    
    % Store results for this task
    allAccuracies.(selName) = struct( ...
        'LDA', nanmean(acc.LDA), ...
        'Logistic', nanmean(acc.Logistic), ...
        'SVM', nanmean(acc.SVM));
    
    allSelectedFeatures.(selName) = struct( ...
        'indices', finalIdx, ...
        'names', {finalNames}, ...
        'accuracies', finalAcc);

    % *** NEW: Save individual per-task feature files needed by analysis scripts ***
    % This is the critical addition that was missing!
    finalNames = cellstr(finalNames(:)); % ensure cell array of char row vectors
    finalAcc = finalAcc(:);              % ensure column vector
    save(fullfile(featDir, sprintf('final_uncorrelated_features_%s.mat', selName)), 'finalNames', 'finalAcc');

    fprintf('  Results: LDA=%.2f%%, Logistic=%.2f%%, SVM=%.2f%%\n', ...
        allAccuracies.(selName).LDA*100, ...
        allAccuracies.(selName).Logistic*100, ...
        allAccuracies.(selName).SVM*100);
end

%% 8) Generate final bar plot
fprintf('\n=== Generating Final Results ===\n');

% Prepare data for plotting
processedTasks = fieldnames(allAccuracies);
nProcessed = numel(processedTasks);

if nProcessed == 0
    error('No tasks were successfully processed.');
end

% Extract accuracies
LDA_acc = zeros(nProcessed, 1);
Logistic_acc = zeros(nProcessed, 1);
SVM_acc = zeros(nProcessed, 1);

for i = 1:nProcessed
    task = processedTasks{i};
    LDA_acc(i) = allAccuracies.(task).LDA;
    Logistic_acc(i) = allAccuracies.(task).Logistic;
    SVM_acc(i) = allAccuracies.(task).SVM;
end

% Create bar plot
figure('Name','Model Accuracies Across Tasks','Position',[100 100 1200 600]);
X = 1:nProcessed;
width = 0.25;

bar(X - width, LDA_acc, width, 'FaceColor', [0.2 0.4 0.6], 'DisplayName', 'LDA');
hold on;
bar(X, Logistic_acc, width, 'FaceColor', [0.8 0.4 0.2], 'DisplayName', 'Logistic');
bar(X + width, SVM_acc, width, 'FaceColor', [0.9 0.7 0.1], 'DisplayName', 'SVM');

xlabel('Tasks');
ylabel('Accuracy');
title('Model Accuracies Across Tasks');
legend('Location', 'best');
xticks(X);
xticklabels(processedTasks);
xtickangle(45);
grid on;
ylim([0 1]);

%% 9) Generate feature selection table
fprintf('\nSelected Features by Task:\n');
fprintf('========================\n');

% Create comprehensive table
featureTable = {};
rowIdx = 1;

for i = 1:nProcessed
    task = processedTasks{i};
    features = allSelectedFeatures.(task);
    
    for j = 1:numel(features.names)
        featureTable{rowIdx, 1} = task;
        featureTable{rowIdx, 2} = j;
        featureTable{rowIdx, 3} = features.names{j};
        featureTable{rowIdx, 4} = features.accuracies(j);
        featureTable{rowIdx, 5} = features.indices(j);
        rowIdx = rowIdx + 1;
    end
end

% Display table
fprintf('%-20s %-4s %-40s %-8s %-8s\n', 'Task', 'Rank', 'Feature Name', 'Accuracy', 'Index');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:size(featureTable, 1)
    fprintf('%-20s %-4d %-40s %-8.2f %-8d\n', ...
        featureTable{i,1}, featureTable{i,2}, featureTable{i,3}, ...
        featureTable{i,4}, featureTable{i,5});
end

%% 10) Save results
resultsFile = fullfile(featDir, 'master_pipeline_results.mat');
save(resultsFile, 'allAccuracies', 'allSelectedFeatures', 'processedTasks', ...
     'LDA_acc', 'Logistic_acc', 'SVM_acc', 'featureTable', '-v7.3');

% *** NEW: Save the main results file needed by analysis scripts ***  
% This is essential for the analysis scripts to work!
save(fullfile(featDir, 'all_model_accuracies.mat'), 'allAccuracies');

fprintf('\nPipeline complete! Results saved to: %s\n', resultsFile);
fprintf('Processed %d tasks successfully.\n', nProcessed);

% *** NEW: Display summary of generated files ***
fprintf('\n=== Files Generated for Analysis Scripts ===\n');
fprintf('Main accuracy file: all_model_accuracies.mat\n');
fprintf('Individual feature files:\n');
for i = 1:nProcessed
    fprintf('  - final_uncorrelated_features_%s.mat\n', processedTasks{i});
end
fprintf('Analysis scripts can now be run!\n');

%% Helper function (from performance_feature_extraction.m)
function [F, L, meta, names] = process_session(session, Fs, ep_len, filters, baseDir)
    % Initialize outputs
    F = [];
    L = [];
    meta = {};
    names = {};
    cnt = 0;
    firstNames = true;

    for rec = 1:2
        fname = sprintf('S04_%s_wk%d_rs%d.mat', session.day, str2double(session.week(end)), rec);
        fpath = fullfile(baseDir, session.week, session.day, fname);
        if ~exist(fpath,'file'), continue; end
        
        try
            D = load(fpath); 
            EEG = D.EEG;
        catch
            continue;
        end

        % Select first 14 channels and resample if needed
        if size(EEG.data, 1) < 14
            continue;
        end
        
        data = EEG.data(1:14,:);
        if EEG.srate ~= Fs
            data = resample(data', Fs, EEG.srate)';
        end

        % Take central 2 minutes
        total = size(data,2);
        window = 2*60*Fs;
        if total < window
            central = data;
        else
            mid = floor(total/2);
            half = floor(window/2);
            central = data(:, mid-half+1 : mid+half);
        end

        % Number of epochs
        nEp = floor(size(central,2) / ep_len);
        if nEp == 0, continue; end
        
        labs = {EEG.chanlocs(1:14).labels};

        for e = 1:nEp
            cnt = cnt + 1;
            seg = central(:, (e-1)*ep_len+1 : e*ep_len);

            % Extract features & names once
            y = arrayfun(@(ch) seg(ch,:)', 1:14, 'UniformOutput', false);
            try
                [vec, nm, ~] = find_JL_ft(y, Fs, labs, [], {}, filters);
                if firstNames
                    names = nm;
                    firstNames = false;
                end

                % Store
                F(cnt, :) = vec;
                L(cnt, 1) = session.is_high_performance;
                meta{cnt, 1} = struct('week',session.week,'day',session.day,'rec',rec,'epoch',e,'score',session.score,'high',session.is_high_performance);
            catch
                cnt = cnt - 1; % Decrement if feature extraction failed
                continue;
            end
        end
    end

    if cnt == 0
        if ~isempty(names)
            F = zeros(0, numel(names));
        else
            F = zeros(0, 1);
        end
        L = zeros(0,1);
        meta = {};
    else
        F = F(1:cnt, :);
        L = L(1:cnt, :);
        meta = meta(1:cnt);
    end
end