function master_pipeline()
    % MASTER_PIPELINE Run the EEG pipeline end-to-end across all 13 tasks.
    
    %% --- ENSURE FUNCTION PATHS ARE AVAILABLE TO PARALLEL WORKERS ---
    % Stop any running parallel pool so workers get the updated path
    pool = gcp('nocreate');
    if ~isempty(pool)
        delete(pool);
    end
    [scriptPath,~,~] = fileparts(mfilename('fullpath'));      % …/code/scripts
    rootDir = fileparts(fileparts(fileparts(scriptPath)));           % …/S04_analysis
    funcDir = fullfile(rootDir, 'code', 'functions', 'multivariate-timeseries-feature-for-EEG-analysis-main');
    % Debug print to verify path
    disp('funcDir:'), disp(funcDir)
    disp('find_JL_ft.m path:'), disp(fullfile(funcDir, 'Feature evaluation functions', 'find_JL_ft.m'))
    addpath(genpath(fullfile(funcDir,'Feature evaluation functions')));
    addpath(genpath(fullfile(funcDir,'HCTSA classification')));
    addpath(genpath(fullfile(funcDir, 'Feature evaluation functions', 'catch22-master', 'wrap_Matlab')));
    % Start pool after all addpath
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool;
    end
    targetFile = fullfile(funcDir, 'Feature evaluation functions', 'find_JL_ft.m');
    if ~exist(targetFile, 'file')
        error('File not found: %s', targetFile);
    end
    addAttachedFiles(pool, {targetFile});

    %% 0) DIRECTORY SETUP
    outFeatDir  = fullfile(rootDir,'code','data','features');
    if ~exist(outFeatDir,'dir'), mkdir(outFeatDir); end
    
    % Base data directory (absolute)
    baseDataDir = '/Users/andreamagri/Desktop/thesis/S04_analysis/data/processed';
    
    %% 1) ADD HELPER-FUNCTION PATHS
    funcDir = fullfile(rootDir,'code','functions', ...
        'multivariate-timeseries-feature-for-EEG-analysis-main');
    addpath(genpath(fullfile(funcDir,'Feature evaluation functions')));
    addpath(genpath(fullfile(funcDir,'HCTSA classification')));
    
    %% 2) DEFINE PERFORMANCE SERIES
    om    = [ 0.3825, -0.7015,  0.3825, -0.0752,  1.7555, -0.8889,  0.8402,  0.7438, -0.9906, -1.4483];
    omd   = [ 0.3753, -1.0148,  0.9546, -0.8063,  0.0741, -0.5514,  0.9546,  1.8350, -1.0148, -0.8063];
    wm    = [ 0.3446, -1.4165,  0.7848, -0.2106, -1.5888,  0.7083, -0.0957,  0.3446, -0.5360,  1.6654];
    wmd   = [ 0.5697, -2.1105,  0.0208,  0.0208,  0.7312,  0.7312, -0.3344,  1.3953, -1.0449,  0.0208];
    fnat  = [-1.5479, -0.6669,  1.1391, -0.6669, -1.2689,  1.2496,  0.5371,  0.9548,  0.0702,  0.1999];
    fnatd = [ 0.8568,  0.0261,  0.1100,  1.1410, -1.3365,  0.8712,  0.2734,  0.1100,  0.0261, -2.0781];
    mst   = [ 1.5571, -0.4104,  0.7878,  0.0546,  0.5566, -1.3808, -0.6293, -0.5994, -1.1534, -1.4483];
    pca1  = [-0.6360,  3.0049, -1.6104,  0.7023, -0.0384, -0.3482, -1.0221, -2.8583,  2.0115,  0.7947];
    
    scores = struct( ...
      'pca1',pca1,'om',om,'omd',omd,'wm',wm,'wmd',wmd, ...
      'fnat',fnat,'fnatd',fnatd,'mst',mst );
    
    % composites
    scores.composite_om        = (om+omd)/2;
    scores.composite_wm        = (wm+wmd)/2;
    scores.composite_fnat      = (fnat+fnatd)/2;
    scores.composite_immediate = (om+wm+fnat)/3;
    scores.composite_delayed   = (omd+wmd+fnatd)/3;
    
    taskNames = { ...
     'pca1','om','omd','wm','wmd', ...
     'fnat','fnatd','mst', ...
     'composite_om','composite_wm','composite_fnat', ...
     'composite_immediate','composite_delayed'};
    nTasks = numel(taskNames);
    
    %% 3) PREALLOCATE RESULTS
    accMat    = nan(nTasks,3);    % [LDA, Logistic, SVM]
    featCells = cell(nTasks,1);    % selected features per task
    
    %% 4) LOOP THROUGH TASKS
    for t = 1:nTasks
      selName = taskNames{t};
      fprintf('\n===== TASK: %s =====\n', selName);
    
      %% STEP 1: FEATURE EXTRACTION
      all_scores = scores.(selName);
      z_hi = 0.25; z_lo = -0.25;
      keep = all_scores>=z_hi | all_scores<=z_lo;
      perf_lbl = all_scores(keep)>=z_hi;
      all_scores = all_scores(keep);
    
      % build session_info
      session_types = {'pre_day','stim_day3','post_day','follow_up1','follow_up3'};
      session_info = {};
      idx_ptr=1; lbl_ptr=1;
      for wk=1:2
        for d=1:numel(session_types)
          if ~keep(idx_ptr)
            idx_ptr=idx_ptr+1; continue;
          end
          session_info{end+1,1} = struct( ...
            'week', sprintf('week_%d',wk), ...
            'day', session_types{d}, ...
            'score', all_scores(lbl_ptr), ...
            'is_high_performance', perf_lbl(lbl_ptr) );
          idx_ptr=idx_ptr+1; lbl_ptr=lbl_ptr+1;
        end
      end
    
      % FIR filters
      Fs=1000; ep_len=5*Fs;
      bands={[1 4],[4 8],[8 12],[13 30]}; order=64;
      filters=cell(1,numel(bands));
      for i=1:numel(bands)
        filters{i}=fir1(order,bands{i}/(Fs/2),'bandpass',hamming(order+1));
      end
    
      % parallel pool
      try
        if isempty(gcp('nocreate')), parpool; end
        use_parallel=true;
      catch
        warning('No parallel toolbox—sequential');
        use_parallel=false;
      end
    
      % process sessions
      N=numel(session_info);
      session_feats=cell(N,1); session_labels=cell(N,1);
      name_list=cell(N,1); session_meta=cell(N,1);
      if use_parallel
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end
        addAttachedFiles(pool, {fullfile(funcDir, 'Feature evaluation functions', 'find_JL_ft.m')});
        D=parallel.pool.DataQueue; afterEach(D,@(~)fprintf('.'));
        parfor i=1:N
          [F,L,mi,nm]=process_session(session_info{i},Fs,ep_len,filters,baseDataDir);
          session_feats{i}=F; session_labels{i}=L; session_meta{i}=mi; name_list{i}=nm;
          send(D,i);
        end
        fprintf('\n');
      else
        for i=1:N
          [F,L,mi,nm]=process_session(session_info{i},Fs,ep_len,filters,baseDataDir);
          session_feats{i}=F; session_labels{i}=L; session_meta{i}=mi; name_list{i}=nm;
          fprintf('.');
        end
        fprintf('\n');
      end
    
      % concatenate & split
      good=~cellfun(@isempty,session_feats);
      allF=vertcat(session_feats{good});
      allLbls=vertcat(session_labels{good});
      ds=name_list{find(good,1)};
      high_feats=allF(allLbls==1,:);
      low_feats= allF(allLbls==0,:);
      rawFile=fullfile(outFeatDir,sprintf('raw_features_%s.mat',selName));
      save(rawFile,'high_feats','low_feats','ds','session_info','selName','-v7.3');
    
      %% STEP 2: FEATURE SELECTION
      R=load(rawFile,'high_feats','low_feats','ds');
      X=[R.high_feats;R.low_feats];
      Y=[ones(size(R.high_feats,1),1);zeros(size(R.low_feats,1),1)];
      numFeat=size(X,2);
      cfnParams=struct('numClasses',2,'cat_labels',Y,...
                       'whatLossFun','classificationError','whatLossUnits','%');
      fparams=struct('numFeatures',numFeat,'ID',(1:numFeat)','description',{R.ds(:)});
      [topIdx,perFeatureAcc]=topFeatureHierClust(X,cfnParams,fparams,10);
      rankFile=fullfile(outFeatDir,sprintf('selected_features_%s.mat',selName));
      save(rankFile,'topIdx','perFeatureAcc','-v7.3');
    
     % ── suppress global clustering errors ─────────────────────────────
        if size(X,2) >= 3
            prevVis = get(0,'DefaultFigureVisible');
            set(0,'DefaultFigureVisible','off');
            try
                correlationHierClust(X, R.ds, perFeatureAcc);
            catch ME
                warning('%s: %s', ME.identifier, ME.message);
            end
            set(0,'DefaultFigureVisible', prevVis);
        end
    
      % greedy uncorrelated
      corrThresh=0.50; accThresh=61; maxKeep=10; selIdx=[];
      for f=topIdx'
        if perFeatureAcc(f)<accThresh, break; end
        if isempty(selIdx) || all(abs(corr(X(:,f),X(:,selIdx)))<corrThresh)
          selIdx(end+1)=f;
        end
        if numel(selIdx)==maxKeep, break; end
      end
      if numel(selIdx) < 5
          warning('Fewer than 5 features passed the thresholds. Filling up to 5 with top features by accuracy.');
          remaining = setdiff(topIdx, selIdx, 'stable');
          nToAdd = 5 - numel(selIdx);
          toAdd = remaining(1:min(nToAdd, numel(remaining)));
          selIdx = [selIdx(:); toAdd(:)]; % ensure both are column vectors
      end
      finalIdx=selIdx; finalNames=R.ds(finalIdx); finalAcc=perFeatureAcc(finalIdx);
      finalFile=fullfile(outFeatDir,sprintf('final_uncorrelated_features_%s.mat',selName));
      save(finalFile,'finalIdx','finalNames','finalAcc','-v7.3');
    
      % suppress final clustering errors
      if numel(finalIdx) >= 3
          prevVis = get(0,'DefaultFigureVisible');
          set(0,'DefaultFigureVisible','off');
          try
              correlationHierClust(X(:, finalIdx), finalNames, finalAcc);
          catch ME
              warning('%s: %s', ME.identifier, ME.message);
          end
          set(0,'DefaultFigureVisible', prevVis);
      end
    
      %% STEP 3: MODEL TRAINING & EVALUATION
      R1=load(rawFile,'high_feats','low_feats'); Xall=[R1.high_feats;R1.low_feats];
      Yall=[ones(size(R1.high_feats,1),1);zeros(size(R1.low_feats,1),1)];
      R2=load(finalFile,'finalIdx'); Xsel=Xall(:,R2.finalIdx);
      kFold=5; rng(1,'twister'); cvp=cvpartition(Yall,'KFold',kFold);
      models={'LDA','Logistic','SVM'};
      for m=models
        acc.(m{1})=nan(kFold,1); AUC.(m{1})=nan(kFold,1);
        precision.(m{1})=nan(kFold,1); recall.(m{1})=nan(kFold,1);
        f1score.(m{1})=nan(kFold,1); confMats.(m{1})=cell(kFold,1);
        ROCXY.(m{1})=cell(kFold,1);
      end
      for fold=1:kFold
        tr=training(cvp,fold); te=test(cvp,fold);
        Xtr=Xsel(tr,:); Ytr=Yall(tr); Xte=Xsel(te,:); Yte=Yall(te);
        mdlL=fitcdiscr(Xtr,Ytr); [pL,sL]=predict(mdlL,Xte);
        acc.LDA(fold)=mean(pL==Yte); CM=confusionmat(Yte,pL);
        confMats.LDA{fold}=CM; [fpr,tpr,~,AUC.LDA(fold)]=perfcurve(Yte,sL(:,2),1);
        ROCXY.LDA{fold}=[fpr,tpr]; TP=CM(2,2); FP=CM(1,2); FN=CM(2,1);
        precision.LDA(fold)=TP/(TP+FP); recall.LDA(fold)=TP/(TP+FN);
        f1score.LDA(fold)=2*(precision.LDA(fold)*recall.LDA(fold))...
                          /(precision.LDA(fold)+recall.LDA(fold));
        mdlG=fitglm(Xtr,Ytr,'Distribution','binomial','Link','logit');
        sg=predict(mdlG,Xte); pg=double(sg>=0.5);
        acc.Logistic(fold)=mean(pg==Yte); CM=confusionmat(Yte,pg);
        confMats.Logistic{fold}=CM; [fpr,tpr,~,AUC.Logistic(fold)]=perfcurve(Yte,sg,1);
        ROCXY.Logistic{fold}=[fpr,tpr]; TP=CM(2,2); FP=CM(1,2); FN=CM(2,1);
        precision.Logistic(fold)=TP/(TP+FP); recall.Logistic(fold)=TP/(TP+FN);
        f1score.Logistic(fold)=2*(precision.Logistic(fold)*recall.Logistic(fold))...
                                /(precision.Logistic(fold)+recall.Logistic(fold));
        mdlS=fitcsvm(Xtr,Ytr,'KernelFunction','linear','Standardize',true);
        [pS,sS]=predict(mdlS,Xte);
        acc.SVM(fold)=mean(pS==Yte); CM=confusionmat(Yte,pS);
        confMats.SVM{fold}=CM; [fpr,tpr,~,AUC.SVM(fold)]=perfcurve(Yte,sS(:,2),1);
        ROCXY.SVM{fold}=[fpr,tpr]; TP=CM(2,2); FP=CM(1,2); FN=CM(2,1);
        precision.SVM(fold)=TP/(TP+FP); recall.SVM(fold)=TP/(TP+FN);
        f1score.SVM(fold)=2*(precision.SVM(fold)*recall.SVM(fold))...
                          /(precision.SVM(fold)+recall.SVM(fold));
      end
      accMat(t,:)=[mean(acc.LDA),mean(acc.Logistic),mean(acc.SVM)];
      featCells{t}=finalNames;
      accFile=fullfile(outFeatDir,'all_model_accuracies.mat');
      if exist(accFile,'file')
        S=load(accFile,'allAccuracies'); allAccuracies=S.allAccuracies;
      else
        allAccuracies=struct();
      end
      allAccuracies.(selName)=struct('LDA',mean(acc.LDA),'Logistic',mean(acc.Logistic),'SVM',mean(acc.SVM));
      save(accFile,'allAccuracies','-v7.3');
    end
    
    %% 5) FINAL PLOT & TABLE
    figure('Name','Model Accuracies Across Tasks','Color','w');
    bar(accMat);
    set(gca,'XTick',1:nTasks,'XTickLabel',taskNames,'XTickLabelRotation',45);
    ylabel('Accuracy'); ylim([0 .85]);
    legend({'LDA','Logistic','SVM'},'Location','northwest');
    grid on;
    
    T = table(taskNames(:), accMat(:,1), accMat(:,2), accMat(:,3), featCells, ...
        'VariableNames',{'Task','LDA_acc','Logistic_acc','SVM_acc','Selected10'});
    disp(T);
    end
    
    %% Subfunction: process_session
    function [F,L,meta,names] = process_session(session, Fs, ep_len, filters, baseDir)
        F=[];L=[];meta={};names={};cnt=0;firstNames=true;
        for rec=1:2
            fname = sprintf('S04_%s_wk%d_rs%d.mat', session.day, str2double(session.week(end)), rec);
            fpath = fullfile(baseDir, session.week, session.day, fname);
            if ~exist(fpath,'file'), continue; end
            D = load(fpath); EEG = D.EEG; data = EEG.data(1:14,:);
            if EEG.srate~=Fs, data = resample(data',Fs,EEG.srate)'; end
            total = size(data,2); window = 2*60*Fs;
            if total>=window
                mid=floor(total/2); half=floor(window/2);
                central=data(:,mid-half+1:mid+half);
            else
                central=data;
            end
            nEp=floor(size(central,2)/ep_len);
            labs={EEG.chanlocs(1:14).labels};
            for e=1:nEp
                cnt=cnt+1;
                seg=central(:,(e-1)*ep_len+1:e*ep_len);
                y=arrayfun(@(ch)seg(ch,:)',1:14,'uni',false);
                [vec,nm,~]=find_JL_ft(y,Fs,labs,[],{},filters);
                if firstNames, names=nm; firstNames=false; end
                F(cnt,:)=vec;
                L(cnt,1)=session.is_high_performance;
                meta{cnt,1}=struct('week',session.week,'day',session.day,'rec',rec,'epoch',e,'score',session.score,'high',session.is_high_performance);
            end
        end
        if cnt>0
            F=F(1:cnt,:); L=L(1:cnt,:); meta=meta(1:cnt);
        end
    end