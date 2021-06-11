function [timeStampCa, WhereImgIs] = timingCaMovies(folder2process)
file2open = ls([folder2process '*trac*.ibw']);
EyeTracker=[];
time = [];
for f = 1:size(file2open,1)
    filePath = [ folder2process file2open(f,:)];
    fileNumber = filePath(strfind(filePath,'00'):strfind(filePath,'00')+2);
    CaFileName = ls([folder2process '*' fileNumber '*sig*']);
    CaFileName = CaFileName(1,:); fileObj = matfile([folder2process CaFileName]);
    CaNframes(f)  =  size(fileObj,'sig',2);
    EyeTrackerTmp = IBWread(filePath);
    timeAll{f} = (0:EyeTrackerTmp.Nsam - 2)' * EyeTrackerTmp.dx + EyeTrackerTmp.x0;
    if f>1
        timeAll{f} = timeAll{f} +  timeAll{f-1}(end)+EyeTrackerTmp.dx ; % one frame gap
    end
    EyeTrackerTmp    = EyeTrackerTmp.y;
    EyeTracker = [EyeTracker EyeTrackerTmp'];
    time = [time timeAll{f}' ];
end

% Initialize some variable
Threshold = 3;
timeStampCa  = zeros (100000,4);
EventCount = 1;
WhereImgIs = zeros (10,5);
WhereImgIsIndex = 1;

% timeStampCa (i,1) = start of UP; timeStampCa (i,2) = Satrt of Down ; timeStampCa (i,3)
% = duration UP state; timeStampCa (i,4) = Duration DOWN state



for  i= 1:(length(EyeTracker)-1)
    
    if ( ( EyeTracker(i) <= Threshold ) && (EyeTracker(i+ 1)  > Threshold ) )
        timeStampCa  (EventCount,1)  = time(i);
        if (EventCount == 1 )
            WhereImgIs (1,1)  =   time(i);
            WhereImgIs (WhereImgIsIndex,3) = EventCount;
        end
        if ((EventCount>=2) )
            timeStampCa  (EventCount,4) = timeStampCa  (EventCount,1) -  timeStampCa  (EventCount -1 ,2);
            if  (timeStampCa  (EventCount,4)> 0.5)
                WhereImgIs (WhereImgIsIndex+1,1) = time(i);
                WhereImgIs (WhereImgIsIndex+1,3) = EventCount;
                WhereImgIs (WhereImgIsIndex,2) = timeStampCa  (EventCount -1 ,2);
                WhereImgIsIndex = WhereImgIsIndex +1;
            end
            if (timeStampCa  (EventCount -1,3)> 0.5)
                WhereImgIs (WhereImgIsIndex+1,1) = time(i);
                WhereImgIs (WhereImgIsIndex+1,3) = EventCount;
                WhereImgIs (WhereImgIsIndex,2) = timeStampCa  (EventCount -1 ,2);
                WhereImgIsIndex = WhereImgIsIndex +1;
            end
        end
    elseif  ( ( EyeTracker(i)  > Threshold ) && (EyeTracker(i+ 1)  <= Threshold ) )
        timeStampCa  (EventCount,2)  = time(i);
        timeStampCa  (EventCount,3) = timeStampCa  (EventCount,2) -  timeStampCa  (EventCount,1);
        EventCount = EventCount +1;
    end
end

timeStampCa (EventCount :end,:) = [];
WhereImgIs (WhereImgIsIndex,2) = timeStampCa (end,2);

WhereImgIs(WhereImgIsIndex + 1 : end, :) = [];
if (WhereImgIs(1,1) == 0)
    WhereImgIs(1,:) = [] ;
end


% Determine if some short recording of less than 1 minute have been
% performed, then remove them (errors during recording time
WhereImgIs(:,4) = WhereImgIs(:,2) - WhereImgIs(:,1);
WhereImgIs((WhereImgIs(:,4)<60),:) = [];
% 
% for i = 1:f
%     if (CaNframes(i) > ((WhereImgIs(i,2)- WhereImgIs(i,1))/(timeStampCa(30,1)- timeStampCa(29,1)))+100)
%         WhereImgIs(i,2) = WhereImgIs(i,1)+CaNframes(i)*(timeStampCa(30,1)- timeStampCa(29,1));
%     end
% end
