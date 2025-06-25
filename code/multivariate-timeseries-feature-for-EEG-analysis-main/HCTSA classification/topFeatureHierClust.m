function [ifeat,testStat] = topFeatureHierClust(feature_space,cfnParams,fparams,numTopFeatures,varargin)
    % topFeatureHierClust: rank and cluster top discriminative features
    % Uses linear SVM accuracy to score features, then clusters top features
    
    numClasses  = cfnParams.numClasses;
    numFeatures = fparams.numFeatures;
    
    % Define loss and test-statistic functions
    fn_loss     = @(yTest,yPred) 100 * mean(yTest == yPred);
    fn_testStat = @(Xtr,Ytr,Xte,Yte) fn_loss(Yte, predict(fitcsvm(Xtr,Ytr,'KernelFunction','linear'), Xte));
    
    % Display header
    fprintf('Computing performance of %u features with %s...\n', numFeatures, 'linear SVM classifier accuracy');
    
    % Time full feature evaluation
    timer   = tic;
    testStat = giveMeStats(feature_space, cfnParams.cat_labels, true);
    
    % Safe elapsed-time print
    elapsed = toc(timer);
    if exist('BF_TheTime','file')
        timeStr = BF_TheTime(elapsed);
    else
        timeStr = sprintf('%.1f seconds', elapsed);
    end
    fprintf(' Done in %s.\n', timeStr);
    
    % Ensure valid stats
    if all(isnan(testStat))
        error('All feature test statistics are NaN. Check data or loss function.');
    end
    
    % Report mean accuracy
    meanStat = mean(testStat, 'omitnan');
    chanceLevel = 100 / numClasses;
    fprintf('Mean accuracy across %u features = %4.2f%%\n', numFeatures, meanStat);
    fprintf('(Random chance = %4.2f%%)\n', chanceLevel);
    
    % Sort features by descending statistic
    [testStatSorted, idxSort] = sort(testStat, 'descend');
    validMask = ~isnan(testStatSorted);
    testStatSorted = testStatSorted(validMask);
    idxSort        = idxSort(validMask);
    
    % List top features
    nTop = min(numTopFeatures, numel(idxSort));
    for i = 1:nTop
        fi = idxSort(i);
        fprintf('[%u] %s -- %4.2f%%\n', fparams.ID(fi), fparams.description{fi}, testStatSorted(i));
    end
    
    % Prepare clustering of top features
    opIndices = idxSort(1:nTop);
    
    % Compute distance matrix
    if exist('BF_pdist','file')
        Dij = BF_pdist(feature_space(:,opIndices)', 'abscorr');
    else
        % fallback to MATLAB pdist (correlation distance)
        Dij = pdist(feature_space(:,opIndices)', 'correlation');
    end
    
    % Build labels for dendrogram
    labels = arrayfun(@(x) sprintf('[%u] %s (%.2f%%)', fparams.ID(x), fparams.description{x}, testStat(x)), opIndices, 'UniformOutput', false);
    
    % Cluster and plot
    if exist('BF_ClusterDown','file')
        BF_ClusterDown(Dij, 'objectLabels', labels, 'whatDistance', 'abscorr');
    else
        Z = linkage(Dij, 'average');
        dendrogram(Z, 'Labels', labels);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ts = giveMeStats(dataMatrix, labels, verbose)
        % giveMeStats: evaluate each feature in isolation
        ts = nan(numFeatures,1);
        if verbose
            fprintf(' Evaluating %u features...\n', numFeatures);
        end
        for k = 1:numFeatures
            try
                ts(k) = fn_testStat(dataMatrix(:,k), labels, dataMatrix(:,k), labels);
            catch
                ts(k) = NaN;
            end
        end
        if verbose
            fprintf(' Completed feature evaluation.\n');
        end
    end
    
    % Return sorted feature indices and original test statistics
    ifeat = idxSort;
    % testStat variable contains per-feature accuracy
    
    end
    