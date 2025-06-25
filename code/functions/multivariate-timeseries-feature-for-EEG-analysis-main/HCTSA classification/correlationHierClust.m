function group_testStat = correlationHierClust(feature_space,ft_dscrp,testStat)
group_testStat = [];
feature_names = ft_dscrp;
num_ft = size(feature_names,1);

group_feature_space = {};
for i = 1:num_ft
    name_now = feature_names{i};
    group_feature_space{i} = feature_space(:,i);
    group_testStat = [group_testStat testStat(i)];
end

%% 
feature_ID = 1:num_ft;
%% compute correlation matrix
Dij = zeros(num_ft,num_ft);
for i = 1:num_ft
    for j = 1:num_ft
        corr_matrix = abs(corr(group_feature_space{1,i},group_feature_space{1,j}));
        Dij(i,j) = 1-mean(corr_matrix(:));
    end
end
%%
distanceMetric = 'abscorr';
op_ind = 1:num_ft;
makeLabel = @(x) sprintf('[%u] %s (%4.2f%s)',feature_ID(x),feature_names{x},...
    group_testStat(x),'');
objectLabels = arrayfun(@(x)makeLabel(x),op_ind,'UniformOutput',0);

clusterThreshold = 0.5;
BF_ClusterDown_new(Dij,'clusterThreshold',clusterThreshold,...
    'whatDistance',distanceMetric,...
    'objectLabels',objectLabels);