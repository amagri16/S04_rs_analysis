%function GiveMeParams 
function params = GiveMeDefaultParams(labels)

% Set the classifier:
params.whatClassifier = 'svm_linear'; % ('svm_linear', 'knn', 'linear', 'fast_linear')

% Number of repeats of cross-validation (reduce variance due to 'lucky splits')
params.numRepeats = 10;

% Number of classes to classify
% (Assume every class is represented in the data):

cat_labels = categorical(labels);
params.cat_labels = cat_labels;
params.classLabels = categories(cat_labels);
numClasses = length(params.classLabels);
params.numClasses = numClasses;

% Cross validation: number of folds (set to 0 for no CV)
params.numFolds = HowManyFolds(cat_labels,numClasses);

%loss units
params.whatLossUnits = '%';

end 

