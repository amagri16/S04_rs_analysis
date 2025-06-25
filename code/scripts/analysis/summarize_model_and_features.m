
% List of tasks (edit if needed)
taskNames = { ...
 'pca1','om','omd','wm','wmd','fnat','fnatd','mst', ...
 'composite_om','composite_wm','composite_fnat', ...
 'composite_immediate','composite_delayed'};

nTasks = numel(taskNames);
LDA_acc = nan(nTasks,1);
Logistic_acc = nan(nTasks,1);
SVM_acc = nan(nTasks,1);
Selected10_str = cell(nTasks,1);

featDir = fullfile('code','data','features');
accFile = fullfile(featDir, 'all_model_accuracies.mat');

% Load all accuracies
S = load(accFile, 'allAccuracies');
allAccuracies = S.allAccuracies;

for t = 1:nTasks
    name = taskNames{t};
    % Load accuracies
    if isfield(allAccuracies, name)
        LDA_acc(t) = allAccuracies.(name).LDA;
        Logistic_acc(t) = allAccuracies.(name).Logistic;
        SVM_acc(t) = allAccuracies.(name).SVM;
    end
    % Load selected features
    finalFile = fullfile(featDir, sprintf('final_uncorrelated_features_%s.mat', name));
    if exist(finalFile, 'file')
        R = load(finalFile, 'finalNames');
        % Convert to comma-separated string for display
        Selected10_str{t} = strjoin(string(R.finalNames), ', ');
    else
        Selected10_str{t} = '';
    end
end

T = table(taskNames(:), LDA_acc, Logistic_acc, SVM_acc, Selected10_str, ...
    'VariableNames', {'Task','LDA_acc','Logistic_acc','SVM_acc','SelectedFeatures'});

disp(T);

% Optionally, export to CSV or Excel
writetable(T, fullfile(featDir, 'model_summary_results.csv'));
% Or for Excel:
% writetable(T, fullfile(featDir, 'model_summary_results.xlsx'));