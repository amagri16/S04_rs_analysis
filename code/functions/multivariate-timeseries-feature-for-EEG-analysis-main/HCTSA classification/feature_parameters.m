function fparams = feature_parameters(feature_subj_all)

fparams.numFeatures = size(feature_subj_all{1}.ft_vec, 2);
fparams.description = feature_subj_all{1}.ft_dscrp;
fparams.ID = (1:size(feature_subj_all{1}.ft_vec, 2))';
end