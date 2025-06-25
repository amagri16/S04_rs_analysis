%% performance_feature_extraction.m (Updated)
% Extract global EEG features and split into high/low performance
% based on user-selected performance series (PCA1 or individual tasks)

clear; clc;

%% 0) Base data directory (absolute)
baseDataDir = '/Users/andreamagri/Desktop/thesis/S04_analysis/data/processed';

%% 1) Add code paths
addpath(fullfile(baseDataDir,'..','..','code','functions', ...
    'multivariate-timeseries-feature-for-EEG-analysis-main','Feature evaluation functions'));

%% 2) Define raw, already z-scored 10-session vectors (1×10 each)
% Task scores
om   = [ 0.3825, -0.7015,  0.3825, -0.0752,  1.7555,  -0.8889,  0.8402,  0.7438, -0.9906, -1.4483]; %done
omd  = [ 0.3753, -1.0148,  0.9546, -0.8063,  0.0741,  -0.5514,  0.9546,  1.8350, -1.0148, -0.8063]; %done
wm   = [ 0.3446, -1.4165,  0.7848, -0.2106, -1.5888,   0.7083, -0.0957,  0.3446, -0.5360,  1.6654]; %done
wmd  = [ 0.5697, -2.1105,  0.0208,  0.0208,  0.7312,   0.7312, -0.3344,  1.3953, -1.0449,  0.0208]; %done
fnat = [-1.5479, -0.6669,  1.1391, -0.6669, -1.2689,   1.2496,  0.5371,  0.9548,  0.0702,  0.1999]; %done
fnatd= [ 0.8568,  0.0261,  0.1100,  1.1410, -1.3365,   0.8712,  0.2734,  0.1100,  0.0261, -2.0781]; %done
mst  = [ 1.5571, -0.4104,  0.7878,  0.0546,  0.5566,  -1.3808, -0.6293, -0.5994, -1.1534, -1.4483]; %done
% PCA1 composite scores
pca1 = [-0.6360,  3.0049, -1.6104,  0.7023, -0.0384,  -0.3482, -1.0221, -2.8583,  2.0115,  0.7947]; %done

%% 3) Pack all series into one struct
scores = struct( ...
    'pca1',           pca1, ...
    'om',             om, ...
    'omd',            omd, ...
    'wm',             wm, ...
    'wmd',            wmd, ...
    'fnat',           fnat, ...
    'fnatd',          fnatd, ...
    'mst',            mst);

%% 4) Compute composites
scores.composite_om   = (scores.om   + scores.omd)   / 2; %done
scores.composite_wm   = (scores.wm   + scores.wmd)   / 2; %done
scores.composite_fnat = (scores.fnat + scores.fnatd) / 2; %done

% composite scores for immediate tasks and delayed tasks
scores.composite_immediate = (scores.om + scores.wm + scores.fnat) / 3; %done
scores.composite_delayed   = (scores.omd + scores.wmd + scores.fnatd) / 3; 

%% 5) Show available series and prompt for selection
fnames = fieldnames(scores);
fprintf('Available performance series:\n');
for i = 1:numel(fnames)
    fprintf('   • %s\n', fnames{i});
end
selName = input('Type the exact name of the series to use (e.g. ''pca1'', ''om'', ''composite_wm''): ','s');
if ~isfield(scores, selName)
    error('Unknown series "%s". Check spelling.', selName);
end

%% 6) Grab chosen vector and create ±0.5-z labels  ─────────────────────────
all_scores = scores.(selName);

% --- label rules ---------------------------------------------------------
z_hi  =  0.5;                       % ≥ +0.5 → high
z_lo  = -0.5;                       % ≤ –0.5 → low

is_high = all_scores >= z_hi;
is_low  = all_scores <= z_lo;
keep    = is_high | is_low;         % logical 10×1 → sessions to analyse

perf_lbl = is_high(keep);           % 1 = high, 0 = low
all_scores = all_scores(keep);      % drop middle-band scores

%% 7) Build session_info (skip middle-band sessions)  ──────────────────────
session_types = {'pre_day','stim_day3','post_day','follow_up1','follow_up3'};
session_info  = {};
idx           = 1;          % 1…10 original session counter
lab_ptr       = 1;          % pointer into perf_lbl (shorter)

for w = 1:2
    for d = 1:numel(session_types)
        if ~keep(idx)          % ignore z in (–0.5, +0.5)
            idx = idx + 1;
            continue
        end
        session_info{end+1,1} = struct( ...               %#ok<*AGROW>
            'week',                sprintf('week_%d',w), ...
            'day',                 session_types{d},     ...
            'score',               all_scores(lab_ptr),  ...
            'is_high_performance', perf_lbl(lab_ptr) );
        idx     = idx + 1;
        lab_ptr = lab_ptr + 1;
    end
end

%% 8) Parameters & FIR filters for PAC
Fs_target    = 1000;
epoch_length = 5 * Fs_target;
bands        = {[1 4],[4 8],[8 12],[13 30]};
order        = 64;
filters      = cell(1,numel(bands));
for i = 1:numel(bands)
    f_nyq = Fs_target/2;
    bn    = bands{i}/f_nyq;
    filters{i} = fir1(order, bn, 'bandpass', hamming(order+1));
end

%% 9) Initialize parallel pool
try
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    use_parallel = true;
catch
    warning('Parallel toolbox unavailable; using sequential mode.');
    use_parallel = false;
end

%% 10) Process all sessions
N              = numel(session_info);
session_feats  = cell(N,1);
session_labels = cell(N,1);
session_meta   = cell(N,1);
name_list      = cell(N,1);

if use_parallel
    D = parallel.pool.DataQueue;
    afterEach(D, @(~) fprintf('.'));
    parfor i = 1:N
        [F,L,mi,nl] = process_session(session_info{i}, Fs_target, epoch_length, filters, baseDataDir);
        session_feats{i}  = F;
        session_labels{i} = L;
        session_meta{i}   = mi;
        name_list{i}      = nl;
        send(D, i);
    end
    fprintf('\n');
else
    for i = 1:N
        [F,L,mi,nl] = process_session(session_info{i}, Fs_target, epoch_length, filters, baseDataDir);
        fprintf('.');
        session_feats{i}  = F;
        session_labels{i} = L;
        session_meta{i}   = mi;
        name_list{i}      = nl;
    end
    fprintf('\n');
end

%% 11) Concatenate & split by label
gthr       = ~cellfun(@isempty, session_feats);
allF       = vertcat(session_feats{gthr});
allLbls    = vertcat(session_labels{gthr});
allMeta    = vertcat(session_meta{gthr});
ds         = name_list{ find(~cellfun(@isempty,name_list),1) };

high_feats = allF(allLbls==1, :);
low_feats  = allF(allLbls==0, :);
sessionsDone = sum(cellfun(@(F) ~isempty(F), session_feats));

datadir = 'code/data/features'; if ~exist(datadir,'dir'), mkdir(datadir); end
save(fullfile(datadir,'raw_features.mat'), ...
     'high_feats','low_feats','ds','allMeta','session_info','selName','-v7.3');

fprintf(['Feature extraction complete:\n' ...
         ' High: %d epochs\n' ...
         ' Low: %d epochs\n' ...
         ' Features: %d\n' ...
         ' Sessions processed: %d/%d\n'], ...
        size(high_feats,1), ...
        size(low_feats,1), ...
        numel(ds), ...
        sessionsDone, N);

%% Helper function
function [F, L, meta, names] = process_session(session, Fs, ep_len, filters, baseDir)
    % Initialize outputs
    F     = [];
    L     = [];
    meta  = {};
    names = {};
    cnt   = 0;
    firstNames = true;

    for rec = 1:2
        fname = sprintf('S04_%s_wk%d_rs%d.mat', session.day, str2double(session.week(end)), rec);
        fpath = fullfile(baseDir, session.week, session.day, fname);
        if ~exist(fpath,'file'), continue; end
        D = load(fpath); EEG = D.EEG;

        % Select first 14 channels and resample if needed
        data = EEG.data(1:14,:);
        if EEG.srate ~= Fs
            data = resample(data', Fs, EEG.srate)';
        end

        % Take central 2 minutes
        total  = size(data,2);
        window = 2*60*Fs;
        if total < window
            central = data;
        else
            mid  = floor(total/2);
            half = floor(window/2);
            central = data(:, mid-half+1 : mid+half);
        end

        % Number of epochs
        nEp  = floor(size(central,2) / ep_len);
        labs = {EEG.chanlocs(1:14).labels};

        for e = 1:nEp
            cnt = cnt + 1;
            seg = central(:, (e-1)*ep_len+1 : e*ep_len);

            % Extract features & names once
            y = arrayfun(@(ch) seg(ch,:)', 1:14, 'UniformOutput', false);
            [vec, nm, ~] = find_JL_ft(y, Fs, labs, [], {}, filters);
            if firstNames
                names = nm;
                firstNames = false;
            end

            % Store
            F(cnt, :)    = vec;
            L(cnt, 1)    = session.is_high_performance;
            meta{cnt, 1} = struct('week',session.week,'day',session.day,'rec',rec,'epoch',e,'score',session.score,'high',session.is_high_performance);
        end
    end

    if cnt == 0
        F    = zeros(0, numel(vec));
        L    = zeros(0,1);
        meta = {};
    else
        F    = F(1:cnt, :);
        L    = L(1:cnt, :);
        meta = meta(1:cnt);
    end
end
