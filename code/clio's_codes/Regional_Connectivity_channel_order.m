function [dscrp, features] = Regional_Connectivity_channel_order(dscrp, features)

% make sure the format is correct
unique_features_single = string();
unique_features_undirected = string();
unique_features_directed = string();

unique_regions = string();

if size(dscrp,1)<size(dscrp,2)
    dscrp = dscrp';
end

if size(features,2)~=size(dscrp,1)
    features = features';
end
%% First extract all single channel feature, pair wise feature
for ft = 1:size(dscrp,1)
    ft_dscrp = dscrp(ft,1);
    ft_dscrp = ft_dscrp{1,1};

    elements = strsplit(ft_dscrp, '_');

    % detect whether pair-wise or single channel
    type = "single";
    for element = elements
        if element == "vs"
             type = "undirected";
        elseif element == "to"
             type = "directed";
        end
    end

    switch type
        case "single"
            if size(elements,2) == 2
                ft_name = elements(2);
            else
                ft_name = strjoin(elements(2:end),'_');
            end

            if ~ismember(ft_name,unique_features_single)

                if unique_features_single == ""
                    unique_features_single = ft_name;
                else
                    unique_features_single = [unique_features_single; ft_name];
                end
            end

            channel_name = regexprep(char(elements(1)), '[^A-Z]', '');
            if ~ismember(channel_name,unique_regions) % looking at the first char, i.e., F4 would be F
                if channel_name == ""
                    unique_regions = channel_name;
                else
                    unique_regions = [unique_regions; channel_name];
                end
            end

        case "undirected"
            if size(elements,2) == 4
                ft_name = elements(4);
            else
                ft_name = strjoin(elements(4:end),'_');
            end

            if ~ismember(ft_name,unique_features_undirected)

                if unique_features_undirected == ""
                    unique_features_undirected = ft_name;
                else
                    unique_features_undirected = [unique_features_undirected; ft_name];
                end
            end

            channel_name = regexprep(char(elements(1)), '[^A-Z]', '');
            if ~ismember(channel_name,unique_regions) % looking at the first char, i.e., F4 would be F
                if channel_name == ""
                    unique_regions = channel_name;
                else
                    unique_regions = [unique_regions; channel_name];
                end
            end

         case "directed"
            if size(elements,2) == 4
                ft_name = elements(4);
            else
                ft_name = strjoin(elements(4:end),'_');
            end

            if ~ismember(ft_name,unique_features_directed)

                if unique_features_directed == ""
                    unique_features_directed = ft_name;
                else
                    unique_features_directed = [unique_features_directed; ft_name];
                end
            end

            channel_name = regexprep(char(elements(1)), '[^A-Z]', '');
            if ~ismember(channel_name,unique_regions) % looking at the first char, i.e., F4 would be F
                if channel_name == ""
                    unique_regions = channel_name;
                else
                    unique_regions = [unique_regions; channel_name];
                end
            end
    end
end
unique_regions(1) = [];

%% Pre-allocate the struct for all combinations

% single channel
ft_name_single = [];
num_unique_region = size(unique_regions,1);

for i = 1:size(unique_regions,1)
    for j = 1:size(unique_features_single,1)
        ft_name_single = [ft_name_single; strcat(unique_regions(i),'_',unique_features_single(j))];
    end
end

ft_all_single = struct();
for i = 1:size(ft_name_single,1)
    ft_all_single.(ft_name_single(i)) = [];
end

% undirected
ft_name_undirected = [];
combinations = nchoosek(1:size(unique_regions,1), 2);
for m = 1:num_unique_region
    combinations = [combinations; m*ones(1,2)];
end

for i = 1:size(combinations,1)
    for j = 1:size(unique_features_undirected, 1)
        ch1 = unique_regions(combinations(i,1));
        ch2 = unique_regions(combinations(i,2));
        ft_name_undirected = [ft_name_undirected; strcat(ch1,'_vs_',ch2,'_',unique_features_undirected(j))];
    end
end

ft_all_undirected = struct();
for i = 1:size(ft_name_undirected,1)
    ft_all_undirected.(ft_name_undirected(i)) = [];
end

% directed
ft_name_directed = [];
for i = 1:size(combinations,1)
    for j = 1:size(unique_features_directed, 1)
        ch1 = unique_regions(combinations(i,1));
        ch2 = unique_regions(combinations(i,2));
        if ch1~=ch2
            ft_name_directed = [ft_name_directed; strcat(ch1,'_to_',ch2,'_',unique_features_directed(j))];
            ft_name_directed = [ft_name_directed; strcat(ch2,'_to_',ch1,'_',unique_features_directed(j))];
        else
            ft_name_directed = [ft_name_directed; strcat(ch1,'_to_',ch2,'_',unique_features_directed(j))];
        end
    end
end
ft_all_directed = struct();
for i = 1:size(ft_name_directed,1)
    ft_all_directed.(ft_name_directed(i)) = [];
end

%%
for ft = 1:size(dscrp,1)

    ft_val = features(:,ft);
    ft_dscrp = dscrp(ft,1);
    ft_dscrp = ft_dscrp{1,1};
    elements = strsplit(ft_dscrp, '_');

    % detect whether pair-wise or single channel
    type = "single";
    for element = elements
        if element == "vs" 
            type = "undirected";
        elseif element == "to"
            type = "directed";
        end
    end

    switch type
        case "single"
            % Extract regional information 
            ch = char(elements(1));
            ch = regexprep(ch, '[^A-Z]', '');
            ft_name = strjoin([ch elements(2:end)],'_');
            ft_mat = ft_all_single.(ft_name);


            % concatonate new feature values
            ft_mat = [ft_mat, ft_val];
            ft_all_single.(ft_name) = ft_mat;


        case "undirected"
            ch1 = char(elements(1));
            ch1 = regexprep(ch1, '[^A-Z]', '');
            ch2 = char(elements(3));
            ch2 = regexprep(ch2, '[^A-Z]', '');

            % ft_name = strjoin([ch1 'vs' ch2 elements(4:end)],'_');
            % ft_mat = ft_all_undirected.(ft_name);
            
            try
                ft_name = strjoin([ch1 'vs' ch2 elements(4:end)],'_');
                ft_mat = ft_all_undirected.(ft_name);
            catch ME
                ft_name = strjoin([ch2 'vs' ch1 elements(4:end)],'_');
                ft_mat = ft_all_undirected.(ft_name);
            end

            % concatonate new feature values
            ft_mat = [ft_mat, ft_val];
            ft_all_undirected.(ft_name) = ft_mat;


        case "directed"       
            ch1 = char(elements(1));
            ch1 = regexprep(ch1, '[^A-Z]', '');
            ch2 = char(elements(3));
            ch2 = regexprep(ch2, '[^A-Z]', '');


            ft_name = strjoin([ch1 'to' ch2 elements(4:end)],'_');
            ft_mat = ft_all_directed.(ft_name);


            % try
            %     ft_name = strjoin([ch1 'to' ch2 elements(4:end)],'_');
            %     ft_mat = ft_all_directed.(ft_name);
            % catch ME
            %     ft_name = strjoin([ch2 'to' ch1 elements(4:end)],'_');
            %     ft_mat = ft_all_directed.(ft_name);
            % end

            % concatonate new feature values
            ft_mat = [ft_mat, ft_val];
            ft_all_directed.(ft_name) = ft_mat;

    end
end

%% Now average and concatonate features into a matrice and their corresponding description

ft_mat_all = [];
ft_dscrp = {};

fields= fieldnames(ft_all_single);

%%############################
% 
% Extract feature types from the field names
featureTypes = cellfun(@(x) strsplit(x, '_'), fields, 'UniformOutput', false);
featureTypes = cellfun(@(x) x{1}, featureTypes, 'UniformOutput', false);

% Create a table to keep track of field names and their feature types
fieldTable = table(fields, featureTypes, 'VariableNames', {'FieldName', 'FeatureType'});

% Sort the table based on the feature type column
sortedFieldTable = sortrows(fieldTable, 'FeatureType');

% Create a new sorted struct
sortedStruct = struct();
for i = 1:height(sortedFieldTable)
    sortedStruct.(sortedFieldTable.FieldName{i}) = ft_all_single.(sortedFieldTable.FieldName{i});
end

field_names = fieldnames(sortedStruct);  

% %%############################


for i = 1:length(field_names)
    field_name = field_names{i};
    ft_mat = sortedStruct.(field_name); 
    avg_ft_mat = mean(ft_mat,2);

    if ~isempty(avg_ft_mat)
        ft_mat_all = [ft_mat_all avg_ft_mat];
        ft_dscrp{end+1,1} = field_name;
    end
end

%%############################

fields = fieldnames(ft_all_undirected);

% Extract feature types from the field names
featureTypes = cellfun(@(x) strsplit(x, '_'), fields, 'UniformOutput', false);

featureTypes_new = {};
for i = 1:length(featureTypes)
    if length(featureTypes{i,1})>=5 && featureTypes{i,1}{5} == "coherence"
        featureTypes_new{i,1}=featureTypes{i,1}{5};

    elseif featureTypes{i,1}{4} == "half"
        featureTypes_new{i,1}=featureTypes{i,1}{6};
    else
        featureTypes_new{i,1}=featureTypes{i,1}{4};
    end
end

featureTypes = featureTypes_new;

%featureTypes = cellfun(@(x) x{4}, featureTypes, 'UniformOutput', false);

% Create a table to keep track of field names and their feature types
fieldTable = table(fields, featureTypes, 'VariableNames', {'FieldName', 'FeatureType'});

% Sort the table based on the feature type column
sortedFieldTable = sortrows(fieldTable, 'FeatureType');

% Create a new sorted struct
sortedStruct = struct();
for i = 1:height(sortedFieldTable)
    sortedStruct.(sortedFieldTable.FieldName{i}) = ft_all_undirected.(sortedFieldTable.FieldName{i});
end

%%############################

field_names = fieldnames(sortedStruct); 
for i = 1:length(field_names)
    field_name = field_names{i};
    ft_mat = sortedStruct.(field_name); 
    avg_ft_mat = mean(ft_mat,2);
    if ~isempty(avg_ft_mat)
        ft_mat_all = [ft_mat_all avg_ft_mat];
        ft_dscrp{end+1,1} = field_name;
    end
end

%%############################

fields = fieldnames(ft_all_directed);

% Extract feature types from the field names
featureTypes = cellfun(@(x) strsplit(x, '_'), fields, 'UniformOutput', false);
featureTypes = cellfun(@(x) x{4}, featureTypes, 'UniformOutput', false);

% Create a table to keep track of field names and their feature types
fieldTable = table(fields, featureTypes, 'VariableNames', {'FieldName', 'FeatureType'});

% Sort the table based on the feature type column
sortedFieldTable = sortrows(fieldTable, 'FeatureType');

% Create a new sorted struct
sortedStruct = struct();
for i = 1:height(sortedFieldTable)
    sortedStruct.(sortedFieldTable.FieldName{i}) = ft_all_directed.(sortedFieldTable.FieldName{i});
end

%%############################

field_names = fieldnames(sortedStruct); 

for i = 1:length(field_names)
    field_name = field_names{i};
    ft_mat = sortedStruct.(field_name); 
    avg_ft_mat = mean(ft_mat,2);
    if ~isempty(avg_ft_mat)
        ft_mat_all = [ft_mat_all avg_ft_mat];
        ft_dscrp{end+1,1} = field_name;
    end
end
dscrp = ft_dscrp;
features = ft_mat_all;