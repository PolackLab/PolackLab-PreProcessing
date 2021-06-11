function [baselineIdx, trialIdx, traceIdx] = getTrialIdx(stimMatrix,timeStamp,trial)
if size(timeStamp,2)>1 && ~iscolumn(timeStamp)
timeStamp=timeStamp(:,1);
end
baselineIdx = find(timeStamp >stimMatrix(trial,1)-1.5 & timeStamp <stimMatrix(trial,1));
baselineIdx = [ baselineIdx(1)-1 ; baselineIdx  ]; % add one frame to be consistent with previous method
trialIdx = find(timeStamp > stimMatrix(trial,1) & timeStamp <stimMatrix(trial,2));
trialIdx = [trialIdx(1)-1 ; trialIdx ]; % add one frame to be consistent with previous method
traceIdx = find(timeStamp >stimMatrix(trial,1)-1.5 & timeStamp <stimMatrix(trial,2)+9);
traceIdx = [traceIdx(1)-1 ; traceIdx ]; % add one frame to be consistent with previous method
end
