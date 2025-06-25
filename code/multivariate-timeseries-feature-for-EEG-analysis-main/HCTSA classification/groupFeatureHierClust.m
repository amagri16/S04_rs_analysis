function group_testStat = groupFeatureHierClust(ft_idx,testStat,feature_space)
group_testStat = [];
feature_names = fieldnames(ft_idx);
num_groups = size(feature_names,1);

group_feature_space = {};
for i = 1:num_groups
    name_now = feature_names{i};
    index = ft_idx.(name_now);
    group_feature_space{i} = feature_space(:,[index(1):index(2)]);
    group_testStat = [group_testStat mean(testStat(index(1):index(2)))];
end

%% 
statUnit = '%';
feature_ID = 1:num_groups;
%% compute correlation matrix
Dij = zeros(num_groups,num_groups);
for i = 1:num_groups
    for j = 1:num_groups
        corr_matrix = abs(corr(group_feature_space{1,i},group_feature_space{1,j}));
        Dij(i,j) = 1-mean(corr_matrix(:));
    end
end
%%
distanceMetric = 'abscorr';
op_ind = 1:num_groups;
makeLabel = @(x) sprintf('[%u] %s (%4.2f%s)',feature_ID(x),feature_names{x},...
    group_testStat(x),statUnit);
objectLabels = arrayfun(@(x)makeLabel(x),op_ind,'UniformOutput',0);

clusterThreshold = 0.5;
BF_ClusterDown_new(Dij,'clusterThreshold',clusterThreshold,...
    'whatDistance',distanceMetric,...
    'objectLabels',objectLabels);