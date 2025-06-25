% List of tasks (edit if needed)
taskNames = { ...
 'pca1','om','omd','wm','wmd','fnat','fnatd','mst', ...
 'composite_om','composite_wm','composite_fnat', ...
 'composite_immediate','composite_delayed'};

featDir = fullfile('code','data','features');

% Gather all selected features and their accuracies
allFeatures = {};
allAccuracies = [];
for t = 1:numel(taskNames)
    finalFile = fullfile(featDir, sprintf('final_uncorrelated_features_%s.mat', taskNames{t}));
    if exist(finalFile, 'file')
        R = load(finalFile, 'finalNames', 'finalAcc');
        nFeat = numel(R.finalNames);
        allFeatures = [allFeatures; string(R.finalNames(:))];
        allAccuracies = [allAccuracies; R.finalAcc(:)];
    end
end

% Unique features and their indices
[uniqueFeatures, ~, idx] = unique(allFeatures);

% Count occurrences and collect accuracies
counts = accumarray(idx, 1);
meanAcc = accumarray(idx, allAccuracies, [], @mean);
maxAcc  = accumarray(idx, allAccuracies, [], @max);
minAcc  = accumarray(idx, allAccuracies, [], @min);

% Sort by frequency (or by meanAcc if you prefer)
[~, sortIdx] = sortrows([counts, meanAcc], [-1 -2]); % sort by count, then meanAcc
featuresSorted = uniqueFeatures(sortIdx);
countsSorted = counts(sortIdx);
meanAccSorted = meanAcc(sortIdx);
maxAccSorted = maxAcc(sortIdx);
minAccSorted = minAcc(sortIdx);

% Plot: Bar for frequency, color for mean accuracy
figure('Name','Feature Frequency and Accuracy Across Tasks','Color','w');
bar(countsSorted, 'FaceColor', [0.2 0.6 0.8]);
hold on;
yyaxis right
plot(meanAccSorted, 'o-r', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Mean Accuracy (%)');
yyaxis left
set(gca, 'XTick', 1:numel(featuresSorted), 'XTickLabel', featuresSorted, 'XTickLabelRotation', 45);
ylabel('Number of Tasks Selected');
xlabel('Feature');
title('Global Feature Importance: Frequency (bar) & Mean Accuracy (red line)');
grid on;
hold off;

% Display table in command window
Tfeat = table(featuresSorted, countsSorted, meanAccSorted, maxAccSorted, minAccSorted, ...
    'VariableNames', {'Feature','NumTasksSelected','MeanAcc','MaxAcc','MinAcc'});
disp(Tfeat);

% Optionally, export to CSV
writetable(Tfeat, fullfile(featDir, 'global_feature_ranking.csv'));