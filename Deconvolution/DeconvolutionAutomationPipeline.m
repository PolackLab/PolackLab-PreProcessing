%Deconvolution Automation Script 2/2019

%% Data loading and Global parameters setting
% file2load = '\\po-ds\John Files\AIM I and II Data for analysis (combined)\Naive Group Data\Individual animal raw data\G78_Naive recording\Naive 1\Impt files for individual background\G78_80517_AV to VA to OT_naive day 1_Processed.mat';
% load(file2load)

%%% check if params where set wherever this script was called from
if ~exist('drift','var')
    drift = 0.01;
end
if ~exist('pnl','var')
    pnl = [0.85 -0.006];
end
%%%

%%Search for the real dt
% dt      = 0.0646; %sampling period
allIntervalls = timeStampCa(:,3)+ timeStampCa(:,4);
allIntervalls = rmoutliers(allIntervalls,'quartiles','thresholdfactor',1.9);
dt = mean(allIntervalls);


Nmovies = size(Signals,2);
Ncells  = size(Signals(end).sig,2);

movSizes=[];
for m=1:Nmovies
    movSizes(m)=size(Signals(m).sig,1);
end
movSizes=[0 movSizes];

DeconvolutionResultsPop = cell(Ncells,1);%preallocating for speed

%% Concatenate all movies

sigAll = sigConcat(Signals);


%% Deal with the NaNs

% find'em
nanidx = find(diff(isnan(sigAll))~=0); %boundaries of NaN territories
if isnan(sigAll(1,1))
    nanidx = [0; nanidx];
end
if mod(length(nanidx),2)~=0 %check if NaN idx went well
    warning('There is not enough NaN boundaries, nanidx{movConds} will be emptied and MLspike will prolly crash');
    nanidx= [];
end
movLength = length(sigAll);

nanLengths = diff(nanidx);
nanLengths = nanLengths(1:2:end);


% cut'em
if ~isempty(nanidx)
    for nanArrays = 1:length(nanidx/2)
        sigAll(isnan(sigAll(:,1)),:) = [];
    end
end



%% Loop and MLspike computations
nProcessed = 0;
parfor n = 1:Ncells % loop on N cells
    tic
    try
        %% Select cell
        calcium = [];
        for movConds = 1:Nmovies   %number of conditions INCLUDING TC
            calcium      = [calcium sigAll(:,n)'];% retrieve calcium signal of cell #n
        end
        calcium      = calcium  /100 + 1; % convert values to percentage and put baseline at 1 instead of 0
        spikesALL = runMLspikeOnOneCell(calcium,nanidx,Nmovies,movSizes,dt,drift,pnl);
        
        %         %% Autocalibration of A, tau and sigma on TC data
        %         parEst                 = spk_autocalibration('par');
        %         parEst.dt              = dt;
        %         parEst.amin            = 0.035;   % set limits for A and tau
        %         parEst.amax            = 0.07;
        %         parEst.taumin          = 0.5;
        %         parEst.taumax          = 1.0;
        %         parEst.driftparam      = drift;
        %         parEst.pnonlin         = pnl ;
        %         % these are the coeffs for the polynomial model used to account for
        %         % Gcamp6f supralinearity; it is the average best values for the
        %         % cells recorded with Gcamp6f
        %
        %         %turn off all displays
        %         parEst.mlspikepar.dographsummary = false;
        %         parEst.display = 'no';
        %
        %         % perform autocal.
        %         [tauest, aest, sigmaest] = spk_autocalibration(calcium,parEst);
        %         %% Param settings after autocal.
        %         par                     = tps_mlspikes('par');
        %         par.dt                  = dt;
        %         par.a                   = aest;
        %         par.tau                 = tauest;
        %         par.finetune.sigma      = sigmaest;
        %         par.algo.nspikemax      = 4;
        %         par.pnonlin             = pnl;
        %         par.drift.parameter     = drift;
        %         par.drift.baselinestart = 1;
        %         par.dographsummary      = false;
        %         par.display             = 'no';
        %
        %
        %
        %
        %         %% MLspike does its magic
        %         spikesALL = cell(1,7);
        %             [spikest, fit, drifting]   = spk_est(calcium,par);
        %             spikesALL{1} = calcium;
        %             spikesALL{2} = spikest;
        %             spikesALL{3} = fit ;
        %             spikesALL{4} = drifting;
        %             spikesALL{5} = [];
        %             spikesALL{6} = par;
        %             spikesALL{7} = length(spikest);
        %
        %
        %
        %         spikesALL = cell2table(spikesALL,'VariableNames',{'Calcium','Spikes','Fit','Drift','Correlation','Parameters','nSpk'});
        %         clear calcium spikest fit drifting
        %
        %
        %         %% Get correlation coeff between real calcium and model
        %
        %             [R,P,RL,RU]        = corrcoef(spikesALL.Calcium, spikesALL.Fit{1});
        %             spikesALL.Correlation{1} = [R(2,1),P(2,1)];
        %
        %
        %         %% Put the NaNs back, adjust spike timing
        %             valuesAdded = 0;
        %             for nanArrays = 1:length(nanidx)/2
        %                 cutPosition = nanidx((nanArrays*2)-1);
        %                 if cutPosition == 1 % if starts with NaNs
        %                     spikesALL.Fit{1} = [NaN(nanLengths(nanArrays),1); spikesALL.Fit{1}];
        %                     spikesALL.Calcium = [NaN(nanLengths(nanArrays),1); spikesALL.Calcium];
        %                     spikesALL.Drift{1} = [NaN(nanLengths(nanArrays),1); spikesALL.Fit{1}];
        %                     spikesALL.Spikes = spikesALL.Spikes + nanLengths(nanArrays)*dt;
        %                     valuesAdded = valuesAdded + nanLengths(nanArrays);
        %                 elseif cutPosition+1+nanLengths(nanArrays) == movLengths % if ends with NaNs
        %
        %                     spikesALL.Fit{1} = [ spikesALL.Fit{1} ; NaN(nanLengths(nanArrays),1)];
        %                     spikesALL.Calcium = [ spikesALL.Calcium ; NaN(nanLengths(nanArrays),1)];
        %                     spikesALL.Drift{1} = [ spikesALL.Drift{1} ; NaN(nanLengths(nanArrays),1)];
        %                     valuesAdded = valuesAdded + nanLengths(nanArrays);
        %
        %                 else % if NaNs in the middle
        %                     spikesALL.Fit{1} = [spikesALL.Fit{1}(1:cutPosition); ...
        %                         NaN(nanLengths(nanArrays),1); ...
        %                         spikesALL.Fit{1}(cutPosition+1:end)];
        %                     spikesALL.Calcium = [spikesALL.Calcium(1:cutPosition); ...
        %                         NaN(nanLengths(nanArrays),1); ...
        %                         spikesALL.Calcium(cutPosition+1:end)];
        %                     spikesALL.Drift{1} = [spikesALL.Drift{1}(1:cutPosition); ...
        %                         NaN(nanLengths(nanArrays),1); ...
        %                         spikesALL.Drift{1}(cutPosition+1:end)];
        %                     spikesALL.Spikes(spikesALL.Spikes>cutPosition*dt) = ...
        %                         spikesALL.Spikes(spikesALL.Spikes>cutPosition*dt) ...
        %                         + nanLengths(nanArrays)*dt;
        %
        %                     valuesAdded = valuesAdded + nanLengths(nanArrays);
        %
        %                 end
        %
        %             end
        %         %% transform the spike times into frames again
        %         spikestmp = spikesALL.Spikes/dt;
        %         [counts, ~ ]= histcounts(spikestmp,0:length(sigAll));
        %         for m = 1:Nmovies
        %             spikesSeparated{m} = counts(sum(movSizes(1:m))+1: sum(movSizes(1:m+1))) ;
        %         end
        %         spikesALL.Spikes = spikesSeparated ;
        %
        %
        %
        %
        %% Saving & cleaning
        DeconvolutionResultsPop{n} = spikesALL;
        %         clear spikesALL
        nProcessed=nProcessed+1;
        fprintf('Deconvolution : %d / %d cells processed\n', n ,size(Signals(end).sig,2))
        
    catch ME
        switch ME.identifier
            case 'MATLAB:dimagree'
                % I would have liked to catch the warning 'no isolated event
                % found!' to stop before but the custom warning  doesnt have
                % a specific Identifier...
                DeconvolutionResultsPop{n} = 'No spike detected';
        end
    end
    toc
end

%% Save .mat file

file2save = sprintf('%s_Deconvolved_drift%1.3f_PNL%.3f%.3f_%s.mat', experimentCode,par.drift.parameter,par.pnonlin,date);
% check if existing, so as not to overwrite
ok=0;
while ~ok
    if isfile(file2save)
        file2save = [file2save(1:end-4) '_bis.mat'];
    else
        ok=1;
    end
end
save(file2save, 'DeconvolutionResultsPop', 'par' )
