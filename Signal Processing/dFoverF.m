function [Signals,Zscore]= dFoverF(Signals,Signals_BG,nMovies,StimMatrix,WhereImgIs,timeStampCa)
% BG fluctuation substraction
for i = 1:nMovies 
    Signals(i).sig = Signals(i).sig - Signals_BG(i).sig + nanmedian(Signals_BG(i).sig);
end
% Create for each neuron a trace including only intertrial intervals
BaselineF = Signals;
for i = 1:size(StimMatrix,1)
    for k=1:nMovies
        if ((StimMatrix(i,1) < WhereImgIs (k,2)) && (StimMatrix(i,2) > WhereImgIs (k,1)))
            MovieNum = k; % Find which movie correspond to the stimulus
        end
    end
    trialIdx = find(timeStampCa (:,1)>StimMatrix(i,1) & timeStampCa (:,1)<StimMatrix(i,2));
    trialIdx = [trialIdx(1)-1 ; trialIdx ]; %consistent with the old one
    BaselineF(MovieNum).sig(trialIdx - WhereImgIs (MovieNum,3), :)= NaN;
end

% Concatenate all the intertrial intervals of all the calcium imaging
% movies
Baseline = [];
for i = 1:nMovies
    Baseline = cat(1,Baseline ,BaselineF(i).sig);
end

% The F0 is determine as the median of the intertrial interval activity.
BaselineF = median(Baseline,1,'omitnan');

% Convert signals in dF/F
for i = 1:nMovies
    Signals(i).sig_raw = Signals(i).sig; %%%%%%%%%%%%%%%%%%%%%%%%%% RAW-BG!!!!!!
    Signals(i).sig_F = BaselineF;
    Signals(i).sig = ((Signals(i).sig_raw - BaselineF )./ BaselineF)*100; % i find this more readable - julien
    
end
BaselinedFoF = ((Baseline- BaselineF) ./ BaselineF)*100; % i find this more readable - julien



% We define z-score as the sdt of the baseline activity
Zscore = std(BaselinedFoF,0,1,'omitnan')';
clear BaselineF BaselinedFoF Baseline c f i k N MovieNum