%% select_uncorrelated_top10.m
% Aggregate per-channel EEG features into global features, rank them via SVM,
% and greedily select the top uncorrelated features.

clear; clc;

%% 1) Setup paths to feature-evaluation functions
[scriptPath,~,~] = fileparts(mfilename('fullpath'));
rootDir   = fullfile(scriptPath,'..','..');
fEvalDir  = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','Feature evaluation functions');
fHCTSADir = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','HCTSA classification');
if exist(fEvalDir,'dir'),    addpath(genpath(fEvalDir));    end
if exist(fHCTSADir,'dir'),   addpath(genpath(fHCTSADir));   end

%% 2) Load raw per-channel features
rawMat = fullfile(rootDir,'code','data','features','raw_features.mat');
if ~exist(rawMat,'file'), error('Cannot find %s', rawMat); end
R = load(rawMat,'high_feats','low_feats','ds');
hf = R.high_feats;
lf = R.low_feats;
ds = R.ds(:);

% Combine into one data matrix and label vector
X = [hf; lf];
Y = [ones(size(hf,1),1); zeros(size(lf,1),1)];

%% 3) Aggregate to global features by averaging across channels
% Extract suffix (text after underscore) for each feature name
tokens   = regexp(ds, '_(.+)$', 'tokens', 'once');
suffixes = vertcat(tokens{:});  % e.g. {'delta'; 'theta';...}
groups   = unique(suffixes);

N = size(X,1);
G = numel(groups);
X_global = zeros(N, G);
for i = 1:G
    cols = find(strcmp(suffixes, groups{i}));
    X_global(:,i) = mean(X(:,cols), 2);
end

ds_global = groups(:);

%% 4) Configure SVM-ranking parameters
numFeat = size(X_global,2);
cfnParams.numClasses    = 2;
cfnParams.cat_labels    = Y;
cfnParams.whatLossFun   = 'classificationError';
cfnParams.whatLossUnits = '%';

fparams.numFeatures = numFeat;
fparams.ID          = (1:numFeat)';
fparams.description = ds_global;

topN = numFeat;

%% 5) Rank all global features with linear SVM
fprintf('Ranking %d global features with linear SVM...\n', numFeat);
[topIdx, perFeatureAcc] = topFeatureHierClust(X_global, cfnParams, fparams, topN);

% Save ranking for record
outRank = fullfile(rootDir,'code','data','features','selected_global_features.mat');
save(outRank, 'topIdx', 'perFeatureAcc', 'ds_global', '-v7.3');
fprintf('Saved global SVM ranking to %s\n', outRank);

%% 6) Plot global correlation clustering
figure;
correlationHierClust(X_global, ds_global, perFeatureAcc);
title('Global feature correlation clustering');

%% 7) Greedy uncorrelated feature selection
corrThresh = 0.5;   % max allowed |r|
accThresh  = 61;    % min accuracy percent
maxKeep    = 10;    % stop after this many

fprintf('\nGreedy selecting global features with acc ≥ %.1f%% and |r| < %.2f (max %d)\n', ...
        accThresh, corrThresh, maxKeep);

selectedIdx = [];
for ii = 1:numFeat
    f   = topIdx(ii);
    acc = perFeatureAcc(f);
    if acc < accThresh
        break;  % accuracy too low
    end
    if isempty(selectedIdx) || all(abs(corr(X_global(:,f), X_global(:,selectedIdx))) < corrThresh)
        selectedIdx(end+1) = f; %#ok<AGROW>
    end
    if numel(selectedIdx) == maxKeep
        break;  % reached desired count
    end
end

finalIdx   = selectedIdx;
finalNames = ds_global(finalIdx);
finalAcc   = perFeatureAcc(finalIdx);

%% 8) Report & save final global features
fprintf('\nSelected %d final global features:\n', numel(finalIdx));
for i = 1:numel(finalIdx)
    fprintf('  [%3d] %-15s — %5.2f%%\n', finalIdx(i), finalNames{i}, finalAcc(i));
end

outFinal = fullfile(rootDir,'code','data','features','final_uncorrelated_global_features.mat');
save(outFinal, 'finalIdx', 'finalNames', 'finalAcc', '-v7.3');
fprintf('Saved final global features to %s\n', outFinal);

%% 9) Plot dependencies among selected global features
figure;
correlationHierClust(X_global(:,finalIdx), finalNames, finalAcc);
title('Dependencies among final selected global features');

fprintf('Global selection pipeline complete.\n');