function [wid,OSI,DSI] = tcProps(tc,fitType,MaxLoc)
% extracts the tuning width, OSI and DSI from the tuning curve
wid = tc.bestPars.wid;

if (contains(fitType,'constant'))
    TCcentered = circshift(tc.bestY, 181 - 0,2);
else
    TCcentered = circshift(tc.bestY, 181 - MaxLoc,2);
end
%OSI = (R pref - R ortho)/(R pref + R ortho)
OSI = (TCcentered(181) - TCcentered(91)) / (abs(TCcentered(181)) + abs(TCcentered(91)));
%DSI = (Rpref - Ropposite)/(Rpref + Ropposite)
DSI = (TCcentered(181) - TCcentered(1)) / (abs(TCcentered(181)) + abs(TCcentered(1)));
