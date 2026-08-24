function samplingAnalysis = tsCopulaSampleJointPeaksMultiVariatePruning_fast( ...
    inputtimestamps,inputtimeseries,varargin)
% Fast, interface-compatible candidate implementation for validation.
% It preserves the original output field names and the symmetric
% strongest-first temporal-envelope pruning concept. The original sampler
% should remain on the MATLAB path while this version is benchmarked.

args.samplingThresholdPrct = [99,99];
args.minPeakDistanceInDaysMonovarSampling = [3,3];
args.maxPeakDistanceInDaysMultivarSampling = 3;
args.marginalAnalysis = cell(1,size(inputtimeseries,2));
args.samplingOrder = 0;
args.peakType = 'allExceedThreshold';
args = tsEasyParseNamedArgs(varargin,args);

t = inputtimestamps(:);
x = inputtimeseries;
nVar = size(x,2);
if size(x,1) ~= numel(t)
    error('Timestamp and time-series lengths differ.');
end

thresholdPrct = expandVector(args.samplingThresholdPrct,nVar, ...
    'samplingThresholdPrct');
minPeakDays = expandVector( ...
    args.minPeakDistanceInDaysMonovarSampling,nVar, ...
    'minPeakDistanceInDaysMonovarSampling');
pairList = nchoosek(1:nVar,2);
nPairs = size(pairList,1);
pairLag = expandVector(args.maxPeakDistanceInDaysMultivarSampling, ...
    nPairs,'maxPeakDistanceInDaysMultivarSampling');
lagMatrix = zeros(nVar);
for ip = 1:nPairs
    lagMatrix(pairList(ip,1),pairList(ip,2)) = pairLag(ip);
    lagMatrix(pairList(ip,2),pairList(ip,1)) = pairLag(ip);
end

dt = tsEvaGetTimeStep(t);
minPeakSamples = minPeakDays/dt;
thresholds = nan(1,nVar);
peakValues = cell(1,nVar);
peakIndices = cell(1,nVar);
peakTimes = cell(1,nVar);

for iv = 1:nVar
    thresholds(iv) = prctile(x(:,iv),thresholdPrct(iv));
    [peakValues{iv},peakIndices{iv}] = findpeaks(x(:,iv), ...
        'minpeakdistance',minPeakSamples(iv));
    peakValues{iv} = peakValues{iv}(:);
    peakIndices{iv} = peakIndices{iv}(:);
    peakTimes{iv} = t(peakIndices{iv});
end

% A global chronological list is used only to choose a unique earliest
% anchor. Every peak is still processed; no peak is consumed at this stage.
allTime = vertcat(peakTimes{:});
allVarCell = arrayfun(@(iv) repmat(iv,numel(peakTimes{iv}),1), ...
    (1:nVar)','UniformOutput',false);
allLocalCell = arrayfun(@(iv) (1:numel(peakTimes{iv}))', ...
    (1:nVar)','UniformOutput',false);
allVar = vertcat(allVarCell{:});
allLocal = vertcat(allLocalCell{:});
allSeriesIndex = vertcat(peakIndices{:});
sortTable = [allTime,allVar,allSeriesIndex,allLocal];
sortTable = sortrows(sortTable,[1 2 3]);

% Collect one block per productive anchor. Repeatedly appending numeric
% matrices caused MATLAB to copy the complete accumulated catalogue at each
% iteration and accounted for more than 90% of the profiled runtime.
nAnchors = size(sortTable,1);
candidateIndicesCell = cell(nAnchors,1);
candidateTimesCell = cell(nAnchors,1);
candidateValuesCell = cell(nAnchors,1);
candidateClassCell = cell(nAnchors,1);
nCandidateBlocks = 0;

for ia = 1:size(sortTable,1)
    anchorTime = sortTable(ia,1);
    anchorVar = sortTable(ia,2);
    anchorSeriesIndex = sortTable(ia,3);
    anchorLocal = sortTable(ia,4);
    lists = cell(1,nVar);

    for iv = 1:nVar
        if iv == anchorVar
            lists{iv} = anchorLocal;
            continue
        end
        tv = peakTimes{iv};
        idxv = peakIndices{iv};
        laterCanonical = tv>anchorTime | ...
            (tv==anchorTime & (iv>anchorVar | ...
            (iv==anchorVar & idxv>anchorSeriesIndex)));
        lists{iv} = find(laterCanonical & ...
            tv-anchorTime<=lagMatrix(anchorVar,iv));
    end
    if any(cellfun(@isempty,lists)); continue; end

    rows = cartesianRows(lists);
    nRows = size(rows,1);
    times = zeros(nRows,nVar);
    values = zeros(nRows,nVar);
    indices = zeros(nRows,nVar);
    for iv = 1:nVar
        times(:,iv) = peakTimes{iv}(rows(:,iv));
        values(:,iv) = peakValues{iv}(rows(:,iv));
        indices(:,iv) = peakIndices{iv}(rows(:,iv));
    end

    validLag = true(nRows,1);
    for ip = 1:nPairs
        i1 = pairList(ip,1); i2 = pairList(ip,2);
        validLag = validLag & ...
            abs(times(:,i1)-times(:,i2))<=pairLag(ip);
    end
    times = times(validLag,:);
    values = values(validLag,:);
    indices = indices(validLag,:);
    if isempty(indices); continue; end

    exceeds = values>=thresholds;
    allExceeds = all(exceeds,2);
    anyExceeds = any(exceeds,2);

    % Preserve the original per-anchor precedence: if at least one fully
    % exceeding combination exists, partial combinations for this anchor
    % are not added. Otherwise, add combinations with any exceedance.
    if any(allExceeds)
        use = allExceeds;
        classNow = ones(nnz(use),1);
    elseif any(anyExceeds)
        use = anyExceeds;
        classNow = 2*ones(nnz(use),1);
    else
        continue
    end
    nCandidateBlocks = nCandidateBlocks+1;
    candidateIndicesCell{nCandidateBlocks} = indices(use,:);
    candidateTimesCell{nCandidateBlocks} = times(use,:);
    candidateValuesCell{nCandidateBlocks} = values(use,:);
    candidateClassCell{nCandidateBlocks} = classNow;
end

if nCandidateBlocks==0
    candidateIndices = zeros(0,nVar);
    candidateTimes = zeros(0,nVar);
    candidateValues = zeros(0,nVar);
    candidateClass = zeros(0,1);
else
    usedBlocks = 1:nCandidateBlocks;
    candidateIndices = vertcat(candidateIndicesCell{usedBlocks});
    candidateTimes = vertcat(candidateTimesCell{usedBlocks});
    candidateValues = vertcat(candidateValuesCell{usedBlocks});
    candidateClass = vertcat(candidateClassCell{usedBlocks});
end

% Defensive de-duplication. The canonical anchor should already make rows
% unique, but unique protects exact-time ties and unusual input grids.
if ~isempty(candidateIndices)
    [candidateIndices,iu] = unique(candidateIndices,'rows','stable');
    candidateTimes = candidateTimes(iu,:);
    candidateValues = candidateValues(iu,:);
    candidateClass = candidateClass(iu,:);
end

samplingOrder = args.samplingOrder;
if any(samplingOrder)
    if numel(samplingOrder)~=2
        error('samplingOrder must be zero or contain two variable indices.');
    end
    realistic = candidateTimes(:,samplingOrder(2))- ...
        candidateTimes(:,samplingOrder(1))>=0;
    candidateIndices = candidateIndices(realistic,:);
    candidateTimes = candidateTimes(realistic,:);
    candidateValues = candidateValues(realistic,:);
    candidateClass = candidateClass(realistic,:);
end

% Same severity definition as the original code. Stable secondary keys make
% ties deterministic without introducing a scientific pivot variable.
severity = mean(candidateValues,2);
if isempty(severity)
    order = zeros(0,1);
else
    tieKeys = [(-severity),min(candidateTimes,[],2),candidateIndices];
    [~,order] = sortrows(tieKeys,1:size(tieKeys,2));
end
candidateIndices = candidateIndices(order,:);
candidateTimes = candidateTimes(order,:);
candidateValues = candidateValues(order,:);
candidateClass = candidateClass(order,:);

% Strongest-first pruning with exactly the temporal-envelope interpretation:
% once [min index,max index] is accepted, intersecting envelopes are rejected.
occupied = false(numel(t),1);
keep = false(size(candidateClass));
for ic = 1:size(candidateIndices,1)
    firstIndex = min(candidateIndices(ic,:));
    lastIndex = max(candidateIndices(ic,:));
    if ~any(occupied(firstIndex:lastIndex))
        keep(ic) = true;
        occupied(firstIndex:lastIndex) = true;
    end
end
candidateIndices = candidateIndices(keep,:);
candidateTimes = candidateTimes(keep,:);
candidateValues = candidateValues(keep,:);
candidateClass = candidateClass(keep,:);

isJoint = candidateClass==1;
jointIndices = candidateIndices(isJoint,:);
jointTimes = candidateTimes(isJoint,:);
jointValues = candidateValues(isJoint,:);
partialIndices = candidateIndices(~isJoint,:);
partialTimes = candidateTimes(~isJoint,:);
partialValues = candidateValues(~isJoint,:);

samplingAnalysis.jointextremes = cat(3,jointTimes,jointValues);
samplingAnalysis.jointextremes2 = cat(3,partialTimes,partialValues);
samplingAnalysis.thresholdsC = thresholds;
% Preserve the original function's legacy single-event orientation (4x1
% rather than 1x4) for strict drop-in compatibility.
if size(jointIndices,1)==1
    samplingAnalysis.jointExtremeIndices = jointIndices.';
else
    samplingAnalysis.jointExtremeIndices = jointIndices;
end
samplingAnalysis.peakIndicesAll = [jointIndices;partialIndices];

marginalAnalysis = args.marginalAnalysis;
samplingAnalysis.jointExtremesNS = [];
if ~isempty(marginalAnalysis)
    nonStatCell = cellfun(@(z) z{2}.nonStatSeries,marginalAnalysis, ...
        'UniformOutput',false);
    nonStat = [nonStatCell{:}];
    if strcmpi(args.peakType,'allExceedThreshold')
        useIndices = jointIndices;
    elseif strcmpi(args.peakType,'anyExceedThreshold')
        useIndices = [jointIndices;partialIndices];
    else
        useIndices = zeros(0,nVar);
    end
    if ~isempty(useIndices)
        ns = zeros(size(useIndices));
        for iv = 1:nVar
            ns(:,iv) = nonStat(useIndices(:,iv),iv);
        end
        samplingAnalysis.jointExtremesNS = ns;
    end
    trendCell = cellfun(@(z) z{2}.trendSeries,marginalAnalysis, ...
        'UniformOutput',false);
    stdCell = cellfun(@(z) z{2}.stdDevSeries,marginalAnalysis, ...
        'UniformOutput',false);
    samplingAnalysis.thresholdsNonStation = cellfun( ...
        @(tr,sd,thr) sd*thr+tr,trendCell,stdCell,num2cell(thresholds), ...
        'UniformOutput',false);
else
    samplingAnalysis.thresholdsNonStation = [];
end

samplingAnalysis.monovariatePeakValues = peakValues;
samplingAnalysis.monovariatePeakTimes = peakTimes;
samplingAnalysis.monovariatePeakIndices = peakIndices;
samplingAnalysis.monovariateExceedanceMask = cell(1,nVar);
samplingAnalysis.monovariateExceedanceValues = cell(1,nVar);
samplingAnalysis.monovariateExceedanceTimes = cell(1,nVar);
samplingAnalysis.monovariateExceedanceIndices = cell(1,nVar);
for iv = 1:nVar
    mask = peakValues{iv}>=thresholds(iv);
    samplingAnalysis.monovariateExceedanceMask{iv} = mask;
    samplingAnalysis.monovariateExceedanceValues{iv} = peakValues{iv}(mask);
    samplingAnalysis.monovariateExceedanceTimes{iv} = peakTimes{iv}(mask);
    samplingAnalysis.monovariateExceedanceIndices{iv} = peakIndices{iv}(mask);
end

fprintf('%d Compound peak events found (fast candidate)\n',size(jointTimes,1));
end

function y = expandVector(x,n,name)
x = x(:).';
if isscalar(x); y = repmat(x,1,n);
elseif numel(x)==n; y = x;
else; error('%s must be scalar or contain %d values.',name,n);
end
end

function rows = cartesianRows(lists)
n = numel(lists);
rows = zeros(1,0);
for iv = 1:n
    values = lists{iv}(:);
    if isempty(rows)
        rows = values;
    else
        nOld = size(rows,1); nNew = numel(values);
        rows = [repelem(rows,nNew,1),repmat(values,nOld,1)]; %#ok<AGROW>
    end
end
end
