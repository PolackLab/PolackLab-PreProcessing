function spikesALL=runMLspikeOnOneCell(calcium,nanidx,Nmovies,movSizes,dt,drift,pnl)
nFrames = length(calcium);
%% Autocalibration of A, tau and sigma on TC data
        parEst                 = spk_autocalibration('par');
        parEst.dt              = dt;
        parEst.amin            = 0.035;   % set limits for A and tau
        parEst.amax            = 0.07;
        parEst.taumin          = 0.5;
        parEst.taumax          = 1.0;
        parEst.driftparam      = drift;
        parEst.pnonlin         = pnl ; 
        % these are the coeffs for the polynomial model used to account for 
        % Gcamp6f supralinearity; it is the average best values for the  
        % cells recorded with Gcamp6f         
        
        %turn off all displays
        parEst.mlspikepar.dographsummary = false;
        parEst.display = 'no';
        
        % perform autocal.
        [tauest, aest, sigmaest] = spk_autocalibration(calcium,parEst);
        %% Param settings after autocal.
        par                     = tps_mlspikes('par');
        par.dt                  = dt;
        par.a                   = aest;
        par.tau                 = tauest;
        par.finetune.sigma      = sigmaest;
        par.algo.nspikemax      = 4;
        par.pnonlin             = pnl;
        par.drift.parameter     = drift; 
        par.drift.baselinestart = 1; 
        par.dographsummary      = false;
        par.display             = 'no';
        


        
        %% MLspike does its magic
        spikesALL = cell(1,7);
            [spikest, fit, drifting]   = spk_est(calcium,par);
            spikesALL{1} = calcium;
            spikesALL{2} = spikest;
            spikesALL{3} = fit ;
            spikesALL{4} = drifting;
            spikesALL{5} = [];
            spikesALL{6} = par;
            spikesALL{7} = length(spikest);
            

        
        spikesALL = cell2table(spikesALL,'VariableNames',{'Calcium','Spikes','Fit','Drift','Correlation','Parameters','nSpk'});
        clear calcium spikest fit drifting
        
        
        %% Get correlation coeff between real calcium and model
        
            [R,P,RL,RU]        = corrcoef(spikesALL.Calcium, spikesALL.Fit{1});
            spikesALL.Correlation{1} = [R(2,1),P(2,1)];
        
        
        %% Put the NaNs back, adjust spike timing
            valuesAdded = 0;
            for nanArrays = 1:length(nanidx)/2
                cutPosition = nanidx((nanArrays*2)-1);
                if cutPosition == 1 % if starts with NaNs
                    spikesALL.Fit{1} = [NaN(nanLengths(nanArrays),1); spikesALL.Fit{1}];
                    spikesALL.Calcium = [NaN(nanLengths(nanArrays),1); spikesALL.Calcium];
                    spikesALL.Drift{1} = [NaN(nanLengths(nanArrays),1); spikesALL.Fit{1}];
                    spikesALL.Spikes = spikesALL.Spikes + nanLengths(nanArrays)*dt;
                    valuesAdded = valuesAdded + nanLengths(nanArrays);
                elseif cutPosition+1+nanLengths(nanArrays) == movLengths % if ends with NaNs
                    
                    spikesALL.Fit{1} = [ spikesALL.Fit{1} ; NaN(nanLengths(nanArrays),1)];
                    spikesALL.Calcium = [ spikesALL.Calcium ; NaN(nanLengths(nanArrays),1)];
                    spikesALL.Drift{1} = [ spikesALL.Drift{1} ; NaN(nanLengths(nanArrays),1)];
                    valuesAdded = valuesAdded + nanLengths(nanArrays);
                    
                else % if NaNs in the middle
                    spikesALL.Fit{1} = [spikesALL.Fit{1}(1:cutPosition); ...
                        NaN(nanLengths(nanArrays),1); ...
                        spikesALL.Fit{1}(cutPosition+1:end)];
                    spikesALL.Calcium = [spikesALL.Calcium(1:cutPosition); ...
                        NaN(nanLengths(nanArrays),1); ...
                        spikesALL.Calcium(cutPosition+1:end)];
                    spikesALL.Drift{1} = [spikesALL.Drift{1}(1:cutPosition); ...
                        NaN(nanLengths(nanArrays),1); ...
                        spikesALL.Drift{1}(cutPosition+1:end)];
                    spikesALL.Spikes(spikesALL.Spikes>cutPosition*dt) = ...
                        spikesALL.Spikes(spikesALL.Spikes>cutPosition*dt) ...
                        + nanLengths(nanArrays)*dt;
                    
                    valuesAdded = valuesAdded + nanLengths(nanArrays);
                    
                end
                
            end
        %% transform the spike times into frames again
        spikestmp = spikesALL.Spikes/dt;
        [counts, ~ ]= histcounts(spikestmp,0:nFrames);
        for m = 1:Nmovies
            spikesSeparated{m} = counts(sum(movSizes(1:m))+1: sum(movSizes(1:m+1))) ;
        end
        spikesALL.Spikes = spikesSeparated ;

        
       
        
        