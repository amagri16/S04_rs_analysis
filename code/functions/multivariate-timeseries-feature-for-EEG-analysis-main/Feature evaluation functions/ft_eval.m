function sampleObj_coll = ft_eval(sampleObj_coll,txt_feature,ifcardiac)

% This function is used to evaluate the features and turn it into a feature
% vector with corresponding descriptions
%
% Author: Junheng Li

%% Log of code

% 29/06/2021, code started
% 21/10/2021, code contined to finish

%% Loading feature description file

fid = fopen(txt_feature,'rt');
A = textscan(fid,'%s%s');% extract the names of all of teh features we want from a text file (contains two  cells) 

ft_names = A{1}; %extract the names of the features form the first cell of A
num_ft = length(ft_names); %define the number of features we are interested in 

%% Feature evaluations

num_samps = length(sampleObj_coll);        % Total number of samples in the collection

% Enables parallel computation 
parfor i = 1:num_samps %iterating through the number of samples 
    
    samp_now = sampleObj_coll{i}; %extract the ith sample structure
    
    y = samp_now.data_EEG; %extract the EEG signal from the ith sample structure
    Fs = samp_now.Fs; %extract the sampling freqeuncy corresponding to this EEG signal
    
    [ft,~] = sleep_ft_extract(y,Fs);        % Core feature evaluation
    ft = struct2cell(ft); % transform the resulting structure into a cell 
    ft = transpose(cell2mat(ft)); %convert the resulting cell into matrix and transpose it 
    
    if ifcardiac              % Possibly add ECG feature evaluation in the future
        
    end
    
    %Checking if the subject has the right snumber of the extracted
    %features 
    if length(ft) ~= num_ft
        samp_now.ft_quality = 0;
        warning(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' does not have the right feature number'])
        continue
    else
        %check if any extracted features of this subject contain NaN values
        samp_now.ft_vec = ft;
        if sum(isnan(ft))>0
            samp_now.ft_quality = 0;
            disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' have NaN feature value'])
            continue
        end
        samp_now.ft_quality = 1;
        samp_now.ft_dscrp = ft_names;
        %If all is okay,siognifiy it by the ft_quality variable = 1 and
        %print out that ft-eval went well 
        disp(['Sample No.',num2str(i),' from Subject No.',num2str(samp_now.Sbj_id),' has been successfully calculated'])
    end
    
    sampleObj_coll{i} = samp_now; %input the extracted features and corresponding iformation back into teh samples form the subject structure 

end




