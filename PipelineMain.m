%Pipeline Main
%% Parameters to input
folder2process = 'D:\Neural Data\B10\110920\'; % The folder where your data is (one session)
useS2P         = 1; % Use s2p or existing hand made segmentation
whatKindOfTC   = {'Spikes'}; % the kind of data you want to use for tuning curve computation ('spikes','AuC','mean')
nTCinDB        = 2; % how many TC blocks are in your database
behavCond      = 'AD';%What label do you want to have for this experiment in the table

%% Some preparation
if useS2P == 1
    useS2Pfiles = 'S2P';
else
    useS2Pfiles = '';
end

experimentCode = ls([folder2process '*001*ball*' ]); % The session code, animal_date
experimentCode = experimentCode(1:10);

%% Stim matrix
totalTime = tic;
stimMatrixTime= tic;
[stimMatrix, IDs, params, colNo] = generateStimMatrix(folder2process); % get the stim matrix
nTrials = size(stimMatrix,1);
stimMatrix(:,14) = stimMatrix(:,1);
stimMatrix(:,1) = [];
stimMatrixTime= toc(stimMatrixTime);
fprintf(' - stimMatrix done in %5.2f seconds \n',stimMatrixTime)
fprintf(' - %i trials found, with the following block structure: %s \n', nTrials, strjoin(IDs))
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)
%% Segmentation
segTime = tic;
if useS2P==1
sigSegSuite2p(folder2process) % do the segmentation using suite2P
end
segTime =toc(segTime);
fprintf(' - Segmentation done (or retrieved) in %5.2f seconds \n',segTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)
%% Df/f and signal prep
% Load the signal files
dffTime=tic;
[Signals,Signals_BG,nMovies,nCells] = loadSigs(folder2process,useS2Pfiles);

% Get the timing of each frame and the different movies location in time
[timeStampCa, WhereImgIs] = timingCaMovies(folder2process);

% Compute Df/f
[Signals,zScore]= dFoverF(Signals,Signals_BG,nMovies,stimMatrix,WhereImgIs,timeStampCa);
clear Signals_BG

dffTime =toc(dffTime);
fprintf(' - Signal files loaded, dF/F computed in %5.2f seconds \n',segTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)

%% Deconvolution
deconvTime =tic;

% to optimize space usage, we capture the variable space before deconv, and
% we will clean the created unused variables after
varsbefore = who; %get names of current variables

%Params setting and script call
drift = 0.01;
pnl = [0.85 -0.006];
files = dir(folder2process);
files = {files(2:end).name};

matches = regexp(files,regexptranslate ('wildcard',sprintf('*_Deconvolved_drift%1.3f_PNL%.3f%.3f_*',drift,pnl)));
matches = ~cellfun(@isempty,matches);
if sum(matches) == 0 % if the deconv file does not exist already
    DeconvolutionAutomationPipeline
    files = dir(folder2process);
    files = {files(2:end).name};
    matches = regexp(files,regexptranslate ('wildcard',sprintf('*_Deconvolved_drift%1.3f_PNL%.3f%.3f_*',drift,pnl)));
    matches = ~cellfun(@isempty,matches);
end
load(files{find(matches,1)}) % load the file saved in the deconv script


% variables cleaning
varsafter = [];
varsnew = [];
varsafter = who;
varsnew = setdiff(varsafter, varsbefore);
clearvars(varsnew{:}, '-except','DeconvolutionResultsPop', 'par')

% formatting the spikes into a variable that has a similar structure than
% Signals, the df/f variable
for neuron = 1:nCells
    if isa(DeconvolutionResultsPop{neuron}, 'table')
        for MovieNum = 1:nMovies
            Spikes(MovieNum).sig(:,neuron) = DeconvolutionResultsPop{neuron}.Spikes{MovieNum};
        end
        % we also extract some deconv stats for the DB
        deconv{neuron} = [DeconvolutionResultsPop{neuron}.Parameters.a, DeconvolutionResultsPop{neuron}.Parameters.tau ...
            , DeconvolutionResultsPop{neuron}.nSpk , DeconvolutionResultsPop{neuron}.Correlation{1}(1,1) ];
        
    else
        for MovieNum = 1:nMovies
            Spikes(MovieNum).sig(:,neuron) = NaN(size(Signals(MovieNum).sig,1),1);
        end
        deconv{neuron} = NaN;
        
    end
end

deconvTime =toc(deconvTime);
fprintf(' - Deconvolution done (or retrieved) in %5.2f seconds \n',deconvTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)

%% Loc, pupil and licking
locTime=tic;
% Locomotion
[ballAll, timeBall] = sbxballmotionPipeline(folder2process,WhereImgIs);
%Licking
[~,~, ~,~, ~, lick, ~, dx] = getIgorsY(folder2process); % get the lick movies
lick = signalCat('join',lick); % concatenate the movies
timeStampLick = linspace(0,(length(lick)-1)*dx,length(lick))'; %have the time of each frame, from the sampling period dx. We want a column

%%%%%%%%%%%%%%%%%%%%%Eye to be integrated

locTime =toc(locTime);
fprintf(' - Locomotion and Pupil done in %5.2f seconds \n',locTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)

%% Cutting the trials (from the different timeseries)
cuttingTime=tic;

% to make it easier, we start by concatenating the calcium sig into one movie
calc = sigConcat(Signals,WhereImgIs,timeStampCa);
spik = sigConcat(Spikes,WhereImgIs,timeStampCa);
   % preallocation for speed
    [dffBaselineMeans, dffMeans,dffAuCs, nSpikes] = deal(NaN(nCells,nTrials));
    [dffTraces, spikesTraces] = deal(num2cell(NaN(nCells,nTrials)));
    [lickTraces,ballTraces] = deal(num2cell(NaN(1,nTrials)));
for trial = 1:nTrials % for each trial
    %% Calcium & Spikes
    % first we get the indices of the trials in the calcium movie
    [baselineIdxNeurons, trialIdxNeurons, traceIdxNeurons] = getTrialIdx(stimMatrix,timeStampCa,trial);

    for cel = 1:nCells  %then we cut baseline, trace AuC and mean from the movie for each cell
        dffBaselineMeans(cel,trial) = nanmean(calc(baselineIdxNeurons,cel));
        dffMeans(cel,trial)= nanmean(calc(trialIdxNeurons,cel))-dffBaselineMeans(cel,trial); %we substract the mean of the prestim baseline
        dffAuCs(cel,trial) = trapz(calc(trialIdxNeurons,cel)-dffBaselineMeans(cel,trial));
        dffTraces{cel,trial} = getByteStreamFromArray(calc(traceIdxNeurons,cel));
        spikesTraces{cel,trial} = getByteStreamFromArray(spik(traceIdxNeurons,cel));
        nSpikes(cel,trial) = sum(spik(trialIdxNeurons,cel));
    end
    
    %% Licking
    
    [~, ~, traceIdx] = getTrialIdx(stimMatrix,timeStampLick,trial);% get the indices of the trial
    lickTraces{trial} = getByteStreamFromArray(lick(traceIdx));% cut the trace from the movie
    
    
    %% Pupil
    
    %% Loc
    [baselineIdx, trialIdx, traceIdx] = getTrialIdx(stimMatrix,timeBall,trial); %get the trial indices
    ballTraces{trial}= getByteStreamFromArray(ballAll(traceIdx)); % cut the trace from the concatenated movie
    
    
end
cuttingTime =toc(cuttingTime);
fprintf(' - Trial slices cut from all variables in %5.2f seconds \n',cuttingTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)
%% Tuning curves
tuningTime=tic;
% prepare matrices with responses to the OT stims
howManyOts = sum(contains(IDs,'OT')); % How many OT blocks
if howManyOts>0 %if there was an OT block
    TCfile = ls([folder2process '*TuningCurves.*']);%check if we already have saved TCs
    if isempty(TCfile) %if not
        
        blocks = stimMatrix(:,13); % indices of all OT trials
        blocks = diff(blocks); % end of block = 1, beginning of block = -1
        blocks = [0; blocks]; %make it equal in size to the n trials
        blocks= [0 (find(blocks)')-1 size(stimMatrix,1)]; %get the boundaries of the different blocks
        otIdx = {};
        otBlock = 0;
        if ~isempty(blocks)
            for b = 1:length(blocks)-1
                if sum(stimMatrix(blocks(b)+1:blocks(b+1),13)~=0)==0
                    otBlock=otBlock+1;
                    otIdx{otBlock} = blocks(b)+1:blocks(b+1);
                end
            end
        end
        if howManyOts~=otBlock
            error('Something went wrong while searching fot the OT blocks')
        end
        
        for ot = 1:howManyOts
            % Keep only trials for orienation tuning curve
            ResponseMatrix = NaN(length(otIdx{ot}),nCells);
            ResponseMatrixSpikes = ResponseMatrix;
            ResponseMatrixAuC = ResponseMatrix;



            if contains(whatKindOfTc,'Means')
                ResponseMatrix(:,2:end) = dffMeans(:,otIdx{ot})';
                ResponseMatrix = [ stimMatrix(otIdx{ot},4),ResponseMatrix] ;

            elseif contains(whatKindOfTc,'Spikes')
                ResponseMatrixSpikes(:,2:end) = nSpikes(:,otIdx{ot})';
                ResponseMatrixSpikes = [ stimMatrix(otIdx{ot},4),ResponseMatrixSpikes] ;

            elseif contains(whatKindOfTc,'AuC')
                ResponseMatrixAuC(:,2:end) = dffAuC(:,otIdx{ot})';
                ResponseMatrixAuC = [ stimMatrix(otIdx{ot},4),ResponseMatrixAuC] ;
                                
            end
            data.ResponseMatrix = ResponseMatrix;
            data.ResponseMatrixSpikes = ResponseMatrixSpikes;
            data.ResponseMatrixAuC =ResponseMatrixAuC;
            clear  trialIdx baselineIdx   MovieNum k i edges j X c Baseline
            
            TCs = fitORtuningCurvesPipeline(data,'DataType',whatKindOfTc);
            allTC{ot} = TCs;
        end
        allTC = [allTC NaN(1,nTCinDB-length(allTC))]; % fill the difference between the recorded TC and the TCs in the database
        save([folder2process experimentCode '_TuningCurves.mat'], 'allTC')
    else %if the file already exists
        load([folder2process TCfile])
    end
    
else %if there is no OT block, fill allTC with NaN
    allTC    = cell(1,nTCinDB);
    allTC(:) = {NaN};
    disp('No OT block')
end

tuningTime =toc(tuningTime);
fprintf(' - Tuning Curves done (or retrieved) in %5.2f seconds \n',tuningTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)


%% Tables asssembly
tableTime =tic;
for t = 1:nTrials % for each trial
    %Trial table
    Trial(t).Trial = t;
    Trial(t).Experiment = experimentCode;
    Trial(t).Behav_Cond = behavCond;
    switch (stimMatrix(t,9))
        case 0
            Trial(t).Block = 'No task';
        case 1
            Trial(t).Block = 'Auditory';
        case 2
            Trial(t).Block = 'AV';
        case 3
            Trial(t).Block = 'Visual';
        case 4
            Trial(t).Block = 'VA';
        case 5
            Trial(t).Block = 'Orientation Tuning';
    end
    
    Trial(t).Visual_Stim = stimMatrix(t,4);
    Trial(t).Auditory_Stim = stimMatrix(t,5);
    Trial(t).Contrast = stimMatrix(t,12);
    Trial(t).Response = stimMatrix(t,7);
    
    switch (stimMatrix(t,8))
        case 0
            Trial(t).Outcome = 'No task';
        case 1
            Trial(t).Outcome = 'Miss';
        case 2
            Trial(t).Outcome = 'FA';
        case 3
            Trial(t).Outcome = 'CR';
        case 4
            Trial(t).Outcome = 'Hit';
    end
    Trial(t).Lick_Trace = lickTraces{t};
    Trial(t).Lick_dt = dx;
    Trial(t).Locomotion_Trace = ballTraces{t};
    
    if stimMatrix(t,13) == 1
        Trial(t).Inclusion = 'included';
    elseif stimMatrix(t,13) == 0
        Trial(t).Inclusion = 'rejected';
    end
    
    for c = 1:nCells
        if t==1
            %cell Table
            Cell_Struct(c).Cell = c;
            Cell_Struct(c).Experiment = experimentCode;
            Cell_Struct(c).Zscore = zScore(c);
            if ~isa(DeconvolutionResultsPop{c},'table')
                Cell_Struct(c).DeconvA_Tau = [];
                Cell_Struct(c).nSpikeTC = NaN;
                Cell_Struct(c).DeconvCorr = NaN;
            else
                Cell_Struct(c).DeconvA_Tau = getByteStreamFromArray(deconv{c}([1,2]));
                Cell_Struct(c).nSpikeTC = (deconv{c}(3));
                Cell_Struct(c).DeconvCorr = (deconv{c}(4));
            end
            
            % To put the TC results in the table we loop on the number of
            % TC blocks
            for tcb = 1:length(allTC) 
                
                for tcType = ['spikes','mean','AuC'] %we will fill the table for each possible type
                    if isfield(allTC{tcb},['TC' tcType ])
                        TuningCurves = allTC{tcb};
                        Cell_Struct(c).(['Tuning_Curve_' tcType '_' num2str(tcb)]) = (getByteStreamFromArray(TuningCurves.(['TC' tcType])(c).bestY'))';
                        Cell_Struct(c).(['Best_Fit_' tcType '_' num2str(tcb)]) =  TuningCurves.(['TC' tcType])(c).bestStr;
                        try
                            if (contains(TuningCurves.(['TC' tcType])(c).bestStr,'constant'))
                                Cell_Struct(c).(['Pref_Orientation_' tcType '_' num2str(tcb)])=NaN;
                            else
                                [~,MaxLoc] = max(TuningCurves.(['TC' tcType])(c).bestY);
                                Cell_Struct(c).(['Pref_Orientation_' tcType '_' num2str(tcb)]) = TuningCurves.(['TC' tcType])(c).bestX(MaxLoc);
                            end
                        catch
                            Cell_Struct(c).(['Pref_Orientation_' tcType '_' num2str(tcb)])=NaN;
                        end
                        
                        try
                            [wid,OSI,DSI] = tcProps(TuningCurves.(['TC' tcType])(c),...
                                TuningCurves.(['TC' tcType])(c).bestStr,...
                                MaxLoc);
                        catch
                            OSI = NaN;
                            DSI=NaN;
                            wid = NaN;
                        end
                        Cell_Struct(c).(['Tuning_Width_' tcType '_' num2str(tcb)]) = wid;
                        Cell_Struct(c).(['DSI_'          tcType '_' num2str(tcb)]) = DSI;
                        Cell_Struct(c).(['OSI_'          tcType '_' num2str(tcb)]) = OSI;
                    else
                        Cell_Struct(c).(['Tuning_Curve_'     tcType '_' num2str(tcb)]) = NaN;
                        Cell_Struct(c).(['Best_Fit_'         tcType '_' num2str(tcb)]) =  '';
                        Cell_Struct(c).(['Pref_Orientation_' tcType '_' num2str(tcb)]) = NaN;
                        Cell_Struct(c).(['Tuning_Width_'     tcType '_' num2str(tcb)]) = NaN;
                        Cell_Struct(c).(['DSI_'              tcType '_' num2str(tcb)]) = NaN;
                        Cell_Struct(c).(['OSI_'              tcType '_' num2str(tcb)]) = NaN;
                        
                    end
                end
                
            end
            load([folder2process ls([folder2process '*CellsStats*'])])
            Cell_Struct(c).isCell = cellsStats.iscell(c);
            Cell_Struct(c).nPix = cellsStats.cellSize(c);

            
            
        end
        %cell trial table
        CellTrial(c + nCells*(t-1)).Experiment  = experimentCode;
        CellTrial(c + nCells*(t-1)).Trial       = t;
        CellTrial(c + nCells*(t-1)).Cell        = c;
        CellTrial(c + nCells*(t-1)).Trace       = dffTraces{c,t};
        CellTrial(c + nCells*(t-1)).TraceSpikes = spikesTraces{c,t};
        CellTrial(c + nCells*(t-1)).nSpikes     = nSpikes(c,t);
        CellTrial(c + nCells*(t-1)).Mean_dFoF   = dffMeans(c,t);
        CellTrial(c + nCells*(t-1)).AuC_dFoF    = dffAuCs(c,t);
    end
end




%% Transform structures into tables to send to the database

Cell_Trial_Table = struct2table(CellTrial);
Trial_Table = struct2table(Trial);
Cell_Table = struct2table(Cell_Struct);

save([folder2process experimentCode '_Cell_Trial_Table.mat'], 'Cell_Trial_Table');
save([folder2process experimentCode '_Cell_Table.mat'], 'Cell_Table');
save([folder2process experimentCode '_Trial_Table.mat'], 'Trial_Table');

tableTime =toc(tableTime);
fprintf(' - Tables done and saved in %5.2f seconds \n',tableTime)
tt = toc(totalTime);
fprintf(' - Total time : %5.2f seconds \n',tt)


%% Upload to SQL

