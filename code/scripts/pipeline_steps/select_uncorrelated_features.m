%% full_feature_selection.m
% Complete pipeline: rank + greedy uncorrelated selection

clear; clc;

%% 1) Add paths to find topFeatureHierClust, correlationHierClust, etc.
[scriptPath,~,~] = fileparts(mfilename('fullpath'));
rootDir   = fullfile(scriptPath,'..','..');
fEvalDir  = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','Feature evaluation functions');
fHCTSADir = fullfile(rootDir,'code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','HCTSA classification');
if exist(fEvalDir,'dir'),    addpath(genpath(fEvalDir));    end
if exist(fHCTSADir,'dir'),   addpath(genpath(fHCTSADir));   end

%% 2) Load raw feature matrices and descriptions
rawMat = fullfile(rootDir,'code','data','features','raw_features.mat');
if ~exist(rawMat,'file'), error('Cannot find raw_features.mat'); end
R = load(rawMat,'high_feats','low_feats','ds');
hf = R.high_feats;
lf = R.low_feats;
ds = R.ds(:);

% assemble full data matrix X and labels Y
X = [hf; lf];             
Y = [ones(size(hf,1),1); zeros(size(lf,1),1)];
numFeat = size(X,2);
if numel(ds)~=numFeat
    error('ds length (%d) ≠ #features (%d)', numel(ds), numFeat);
end

%% 3) Set up SVM‐ranking parameters
cfnParams.numClasses    = 2;
cfnParams.cat_labels    = Y;
cfnParams.whatLossFun   = 'classificationError';
cfnParams.whatLossUnits = '%';

fparams.numFeatures = numFeat;
fparams.ID          = (1:numFeat)';
fparams.description = ds;

%% 4) Rank all features via single‐feature linear SVM
topN = 10;  % only used in the plotting function
fprintf('Ranking features 1…%d with linear SVM…\n', numFeat);
[topIdx, perFeatureAcc] = topFeatureHierClust(X, cfnParams, fparams, topN);

% save ranking for record
outRank = fullfile(rootDir,'code','data','features','selected_features.mat');
save(outRank,'topIdx','perFeatureAcc','-v7.3');
fprintf('Saved SVM ranking to %s\n', outRank);

%% 5) (Optional) Global correlation clustering of all features
fprintf('Plotting global correlation clustering of all features…\n');
figure;
correlationHierClust(X, ds, perFeatureAcc);
title('Global feature correlation clustering');

%% 6) Greedy uncorrelated + high‐accuracy selection
% Settings:
corrThresh = 0.5;    % maximum allowed abs(corr) between any two kept features
accThresh  = 61;     % only keep features with accuracy ≥ this percent
maxKeep    = 10;     % stop after this many survivors

fprintf('\nGreedy selecting features with acc ≥ %.1f%% and mutual |r| < %.2f (max %d)\n', ...
    accThresh, corrThresh, maxKeep);

selectedIdx = [];
for ii = 1:numFeat
    f = topIdx(ii);             % next best feature
    acc = perFeatureAcc(f);
    if acc < accThresh
        break;                  % no more high‐accuracy features down the list
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
        break;                  % reached desired count
    end
end

% After greedy selection
if numel(selectedIdx) < 5
    warning('Fewer than 5 features passed the thresholds. Filling up to 5 with top features by accuracy.');
    % Find features not already selected
    remaining = setdiff(topIdx, selectedIdx, 'stable');
    % Add up to 5 total
    nToAdd = 5 - numel(selectedIdx);
    toAdd = remaining(1:min(nToAdd, numel(remaining)));
    selectedIdx = [selectedIdx(:); toAdd(:)];
end

finalIdx   = selectedIdx;
finalNames = ds(finalIdx);
finalAcc   = perFeatureAcc(finalIdx);

%% 7) Report & save final set
fprintf('\nSelected %d final features:\n', numel(finalIdx));
for i = 1:numel(finalIdx)
    fprintf('  [%3d] %-15s — %5.2f%%\n', ...
        finalIdx(i), finalNames{i}, finalAcc(i));
end

outFinal = fullfile(rootDir,'code','data','features','final_uncorrelated_features.mat');
save(outFinal,'finalIdx','finalNames','finalAcc','-v7.3');
fprintf('Saved final features to %s\n', outFinal);

%% 8) Plot dependencies among the final features
fprintf('Plotting %d×%d dependency matrix…\n', numel(finalIdx), numel(finalIdx));
figure;
correlationHierClust(X(:,finalIdx), finalNames, finalAcc);
title('Dependencies among final selected features');

fprintf('Pipeline complete.\n');