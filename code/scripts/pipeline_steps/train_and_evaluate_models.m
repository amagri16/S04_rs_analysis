%% train_and_evaluate_models.m
% Train & evaluate LDA, logistic regression, and linear SVM
% with extended metrics: accuracy, AUC, precision, recall, F1, and feature weights.

clear; clc;

%% 1) Paths and load data
[scriptPath,~,~] = fileparts(mfilename('fullpath'));
parentDir = fileparts(scriptPath);        % …/code/scripts
rootDir   = fileparts(parentDir);        % …/S04_analysis
featDir   = fullfile(rootDir,'code','data','features');

% Load raw features and labels
rawFile = fullfile(featDir,'raw_features.mat');
assert(exist(rawFile,'file')==2, 'raw_features.mat not found at %s', rawFile);
R_raw = load(rawFile,'high_feats','low_feats');
X_all = [R_raw.high_feats; R_raw.low_feats];
Y_all = [ones(size(R_raw.high_feats,1),1); zeros(size(R_raw.low_feats,1),1)];

% Load selected feature indices and names
selFile = fullfile(featDir,'final_uncorrelated_features.mat');
assert(exist(selFile,'file')==2, 'final_uncorrelated_features.mat not found at %s', selFile);
R_sel = load(selFile,'finalIdx','finalNames');
X_sel = X_all(:, R_sel.finalIdx);
featNames = R_sel.finalNames;

%% 2) Cross-validation setup
kFold = 5;                      % number of folds
rng(1,'twister');               % for reproducibility
cvp = cvpartition(Y_all,'KFold',kFold);
models = {'LDA','Logistic','SVM'};

% Preallocate storage
for m=1:numel(models)
    acc.(models{m}) = zeros(kFold,1);
    AUC.(models{m}) = zeros(kFold,1);
    precision.(models{m}) = zeros(kFold,1);
    recall.(models{m})    = zeros(kFold,1);
    f1score.(models{m})   = zeros(kFold,1);
    confMats.(models{m})  = cell(kFold,1);
    ROCXY.(models{m})      = cell(kFold,1);  % { [FPR,TPR] }
end

%% 3) Train, evaluate, and compute metrics per fold
for fold = 1:kFold
    tr = training(cvp,fold);
    te = test(cvp,fold);
    Xtr = X_sel(tr,:); Ytr = Y_all(tr);
    Xte = X_sel(te,:); Yte = Y_all(te);

    % ---- LDA ----
    mdlLDA = fitcdiscr(Xtr, Ytr);
    [predLDA, scoreLDA] = predict(mdlLDA, Xte);
    acc.LDA(fold) = mean(predLDA==Yte);
    confMats.LDA{fold} = confusionmat(Yte,predLDA);
    [Xl,Yl,~,AUC.LDA(fold)] = perfcurve(Yte,scoreLDA(:,2),1);
    ROCXY.LDA{fold} = [Xl(:),Yl(:)];
    % precision/recall/F1
    TP = confMats.LDA{fold}(2,2);
    FP = confMats.LDA{fold}(1,2);
    FN = confMats.LDA{fold}(2,1);
    precision.LDA(fold) = TP/(TP+FP);
    recall.LDA(fold)    = TP/(TP+FN);
    f1score.LDA(fold)   = 2*(precision.LDA(fold)*recall.LDA(fold))/(precision.LDA(fold)+recall.LDA(fold));

    % ---- Logistic ----
    mdlLog = fitglm(Xtr, Ytr, 'Distribution','binomial','Link','logit');
    scoresLog = predict(mdlLog, Xte);
    predLog = double(scoresLog>=0.5);
    acc.Logistic(fold) = mean(predLog==Yte);
    confMats.Logistic{fold} = confusionmat(Yte,predLog);
    [Xg,Yg,~,AUC.Logistic(fold)] = perfcurve(Yte,scoresLog,1);
    ROCXY.Logistic{fold} = [Xg(:),Yg(:)];
    TP = confMats.Logistic{fold}(2,2);
    FP = confMats.Logistic{fold}(1,2);
    FN = confMats.Logistic{fold}(2,1);
    precision.Logistic(fold) = TP/(TP+FP);
    recall.Logistic(fold)    = TP/(TP+FN);
    f1score.Logistic(fold)   = 2*(precision.Logistic(fold)*recall.Logistic(fold))/(precision.Logistic(fold)+recall.Logistic(fold));

    % ---- SVM ----
    mdlSVM = fitcsvm(Xtr, Ytr, 'KernelFunction','linear','Standardize',true);
    [predSVM, scoreSVM] = predict(mdlSVM, Xte);
    acc.SVM(fold) = mean(predSVM==Yte);
    confMats.SVM{fold} = confusionmat(Yte,predSVM);
    [Xs,Ys,~,AUC.SVM(fold)] = perfcurve(Yte,scoreSVM(:,2),1);
    ROCXY.SVM{fold} = [Xs(:), Ys(:)];
    TP = confMats.SVM{fold}(2,2);
    FP = confMats.SVM{fold}(1,2);
    FN = confMats.SVM{fold}(2,1);
    precision.SVM(fold) = TP/(TP+FP);
    recall.SVM(fold)    = TP/(TP+FN);
    f1score.SVM(fold)   = 2*(precision.SVM(fold)*recall.SVM(fold))/(precision.SVM(fold)+recall.SVM(fold));
end

%% 4) Aggregate metrics and print summary
fprintf('\n=== Cross-Validated Metrics (mean ± std) ===\n');
for m = 1:numel(models)
    name = models{m};
    fprintf('%s: Accuracy=%.2f%%±%.2f%%, AUC=%.3f±%.3f, Prec=%.3f±%.3f, Rec=%.3f±%.3f, F1=%.3f±%.3f\n', ...
        name, mean(acc.(name))*100, std(acc.(name))*100, ...
        mean(AUC.(name)), std(AUC.(name)), ...
        mean(precision.(name)), std(precision.(name)), ...
        mean(recall.(name)), std(recall.(name)), ...
        mean(f1score.(name)), std(f1score.(name)) );
end

%% 5) Plot average confusion matrices
figure('Name','Avg Confusion Matrices','Position',[50 50 1000 300]);
for i = 1:3
    subplot(1,3,i);
    Csum = sum(cat(3,confMats.(models{i}){:}),3);
    confusionchart(Csum, {'Low','High'});
    title(models{i});
end

%% 6) Plot mean ROC curves
figure('Name','Mean ROC Curves','Position',[100 400 600 500]); hold on;
FPR = linspace(0,1,100);
for m = 1:numel(models)
    % interpolate all ROC curves onto a common FPR grid, handling duplicates
    TPRs = NaN(kFold, numel(FPR));
    for f = 1:kFold
        xy = ROCXY.(models{m}){f};
        % ensure unique FPR points for interp1
        [fpr_u, ia] = unique(xy(:,1));
        tpr_u = xy(ia,2);
        % interpolate (NaN for out-of-bounds)
        TPRs(f,:) = interp1(fpr_u, tpr_u, FPR, 'linear', NaN);
    end
    meanTPR = nanmean(TPRs,1);
    plot(FPR, meanTPR, 'LineWidth', 2);
end
xlabel('False Positive Rate'); ylabel('True Positive Rate');
legend(models,'Location','southeast'); title('Mean ROC'); grid on; hold off;

%% 7) Feature weights on full data on full data
% Logistic regression coefficients
mdlLogFull = fitglm(X_sel, Y_all, 'Distribution','binomial','Link','logit');
wLog = mdlLogFull.Coefficients.Estimate(2:end);
% SVM weights
mdlSVMFull = fitcsvm(X_sel, Y_all, 'KernelFunction','linear','Standardize',true);
wSVM = mdlSVMFull.Beta;

figure('Name','Feature Weights','Position',[200 100 800 400]);
subplot(1,2,1);
bar(wLog); xticks(1:numel(featNames)); xticklabels(featNames); xtickangle(45);
title('Logistic Regression Weights'); ylabel('Coefficient');
subplot(1,2,2);
bar(wSVM); xticks(1:numel(featNames)); xticklabels(featNames); xtickangle(45);
title('Linear SVM Weights');

%% 8) Save expanded evaluation results
outFile = fullfile(featDir,'model_evaluation_extended.mat');
save(outFile, 'acc','AUC','precision','recall','f1score','confMats','ROCXY','wLog','wSVM','featNames','cvp');
fprintf('\nSaved extended evaluation results to %s\n', outFile);

% === Save mean accuracies for this task ===
% Load task name from raw_features.mat
R_raw = load(rawFile, 'selName');
taskName = R_raw.selName;

% Prepare accuracy struct
modelAccuracies = struct('LDA', mean(acc.LDA), ...
                         'Logistic', mean(acc.Logistic), ...
                         'SVM', mean(acc.SVM));

% Load existing results if they exist
accFile = fullfile(featDir, 'all_model_accuracies.mat');
if exist(accFile, 'file')
    S = load(accFile);
    allAccuracies = S.allAccuracies;
else
    allAccuracies = struct();
end

% Store/update for this task
allAccuracies.(taskName) = modelAccuracies;

% Save back
save(accFile, 'allAccuracies', '-v7.3');
