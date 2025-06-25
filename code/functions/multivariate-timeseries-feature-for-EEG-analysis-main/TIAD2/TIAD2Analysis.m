

%% Some kind of high level script 

% Main benefit I think is automated directory handling. 
% No need to think about where files go, just define a directory e.g.
% EEG_analysis, and the folders will be made / found. 
% Also avoids unecessary computation if certain files already exist - means
% you can just add to what already exists I guess. 

% Input: Raw files and a directory

% Output: Folders for preprocessed files, epoched data, images, and eventually features. 

% addpath '/Users/simonwilliamson/Desktop/PhD/EEG analysis/eeglab2024.0';  % Should only need to run these when first starting up matlab
% eeglab;
% addpath '/Users/simonwilliamson/Desktop/PhD/EEG analysis/fieldtrip-20240110';
% ft_defaults

function TIAD2Analysis(data_direc, output_direc)

% Overarching comments here
    
    % preprocess:
    Preproc(data_direc, output_direc, '03', 1000, false, true)   % Preproc.m is the script
    % Preproc('D:\Study Data', '/Users/simonwilliamson/Desktop/PhD/EEG analysis/', '03', 1000, false, true)

    




end 