function sig = sigConcat(Signals,WhereImgIs,timeStampCa)
%small function to concat the signals from different movies 


arguments
    
    Signals (1,:) {isa(Signals,'struct')}
    WhereImgIs (:,5) = [] 
    timeStampCa (:,4) = []
end

sig=[];
for n = 1:size(Signals,2)
    if ~isempty(WhereImgIs)
    try
    foundMovieLength = length(WhereImgIs(n,3) : WhereImgIs(n+1,3)-1);
    catch
        foundMovieLength = length(WhereImgIs(n,3) : length(timeStampCa));
    end
    else
        foundMovieLength=[];
    end
    sig = [sig ;Signals(n).sig ; NaN(foundMovieLength-size(Signals(n).sig,1),size(Signals(n).sig,2))];
end