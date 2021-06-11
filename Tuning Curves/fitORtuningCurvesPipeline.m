function TCs = fitORtuningCurvesPipeline(data,options)
arguments
    data (1,:) {isa(data,'struct')}
    options.DataType (1,:) char {mustBeMember(options.DataType,{'Spikes','Mean','AUC'})} = 'Spikes'
end


if iscell(options.DataType)
    nAsked = length(options.DataType);
else
    nAsked=1;
end

if length(data)~= nAsked
    error(['You asked for ' num2str(nAsked) ' TC types but provided ' num2str(length(data)) ' data types'] )
end

curvesToTest={'constant','circular_gaussian_180','circular_gaussian_360','direction_selective_circular_gaussian'};

TCs = struct;


for i = 1:1: (size (data.ResponseMatrix,2))-1
    tic
    
    if ismember(options.DataType,'Mean')
        %Perform the fit using the mean dFoF recorded during the visual stimulus
        
        B=jdBayesPhysV1('deg',data.ResponseMatrix(:,1),'resp',data.ResponseMatrix(:,i+1),'curvenames',curvesToTest,'unit','dfof','parfor',true);
        
        %Write the results in a structure
        
        TCs.TCmean(i).bestStr = B.bestStr;
        TCs.TCmean(i).bestIdx = B.bestIdx;
        TCs.TCmean(i).bestPars = B.bestPars;
        TCs.TCmean(i).pars = B.pars;
        TCs.TCmean(i).bestX = B.bestX;
        TCs.TCmean(i).bestY = B.bestY;
        
        
    elseif ismember(options.DataType,'AUC')
        %Perform the fit using the AUC dFoF recorded during the visual stimulus
        
        C=jdBayesPhysV1('deg',data.AreaMatrix(:,1),'resp',(data.AreaMatrix(:,i+1))*1e-3,'curvenames',curvesToTest,'unit','dfof','parfor',true);
        
        %Write the results in a structure
        
        TCs.TCarea(i).bestStr = C.bestStr;
        TCs.TCarea(i).bestIdx = C.bestIdx;
        TCs.TCarea(i).bestPars = C.bestPars;
        TCs.TCarea(i).pars = C.pars;
        TCs.TCarea(i).bestX = C.bestX;
        TCs.TCarea(i).bestY = C.bestY;
        
    elseif ismember(options.DataType,'Spikes')
        
        %Perform the fit using the mean spiking rate  recorded during the visual stimulus
        if sum(isnan(data.ResponseMatrixSpikes(:,i+1)))==0
            D=jdBayesPhysV1('deg',data.ResponseMatrixSpikes(:,1),'resp',data.ResponseMatrixSpikes(:,i+1),'curvenames',curvesToTest,'unit','spikerate','parfor',true);
            
            %Write the results in a structure
            
            TCs.TCspikes(i).bestStr = D.bestStr;
            TCs.TCspikes(i).bestIdx = D.bestIdx;
            TCs.TCspikes(i).bestPars = D.bestPars;
            TCs.TCspikes(i).pars = D.pars;
            TCs.TCspikes(i).bestX = D.bestX;
            TCs.TCspikes(i).bestY = D.bestY;
        else
            TCs.TCspikes(i).bestStr = NaN;
            TCs.TCspikes(i).bestIdx = NaN;
            TCs.TCspikes(i).bestPars = NaN;
            TCs.TCspikes(i).pars = NaN;
            TCs.TCspikes(i).bestX = NaN;
            TCs.TCspikes(i).bestY = NaN;
        end
        toc
    end
     
    fprintf(' %1.0f / %3.0f Tuning curves computed ',i,size (data.ResponseMatrix,2)-1)
    
    
end