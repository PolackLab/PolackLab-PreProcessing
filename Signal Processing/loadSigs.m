function [Signals,Signals_BG,Nmovies,nCells] = loadSigs(folder2process,useS2Pfiles)
ListSignals = ls([folder2process '*' useS2Pfiles  '*signals.*']); 
NumSigBG = 0;
NumSig = 0;

for i= 1:size(ListSignals,1)
file2open = [folder2process ListSignals(i,:)];    
    if ( ~contains(file2open,'BG'))
        NumSig = NumSig+1; 
        Signals(NumSig) = load(file2open);
        
    else
        NumSigBG = NumSigBG+1;
        Signals_BG(NumSigBG) = load([file2open]);
    end
end
Nmovies=NumSig;
if NumSigBG~=NumSig
error('The number of sig files is different than the number of BG files')
end
nCells = size(Signals(1).sig,2);