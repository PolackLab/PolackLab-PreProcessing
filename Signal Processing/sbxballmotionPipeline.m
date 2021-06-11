function [ballAll, timeAll] = sbxballmotionPipeline(folder2process,WhereImgIs)
%Function that extract the absolute ball movement and time stamps
% Gives a concatenated, single array


% before extracting the values, we need the time gaps between the sbx
% movies
moviesGaps(1)= WhereImgIs(1,1);
for m = 2:size(WhereImgIs,1)
moviesGaps(m) = WhereImgIs(m,1)-WhereImgIs(m-1,2);
end


%Then we compute the data
bFiles =  ls([folder2process '*_ball.mat']);
ballAll = []; 
timeAll = [0];
for b = 1:size(bFiles,1)
    ballabs=[];
    load([folder2process bFiles(b,:)],'-mat');            %
    
    data = squeeze(data);
    M = double(max(data(:)));
    m = double(min(data(:)));
    
    ball = zeros(1,size(data,3));
    
    % estimate threshold
    
    th = zeros(1,100);
    idx = round(linspace(1,size(data,3),100));
    for(i=1:100)
        z = double(data(:,:,i));
        z = (z-m)/(M-m);
        th(i) = graythresh(z);
    end
    th = median(th);
    
    
    % estimate translation...
    
    z1 = double(data(:,:,1));
    z1 = (z1-m)/(M-m);
    z1 = uint8(im2bw(z1,th));
    
    for n=1:size(data,3)-1
        z2 = double(data(:,:,n+1));
        z2 = (z2-m)/(M-m);
        z2 = uint8(im2bw(z2,th));
        [u,v] = fftalign(z1,z2);
        ball(n) = v+1i*u;
        z1 = z2;
        ballabs(n) = sqrt(v^2+u^2);
    end
    ballAll = [ballAll ballabs 0];
    lastTimeVal = timeAll(end);
    timeAll= [timeAll ; time + lastTimeVal + moviesGaps(b) ];
    
    
end
timeAll(1)=[];

% Some potential problem:
% There's a difference between the calcium frames and the ball frames times that 
% accumulates, but it ends up being less than 4 ms qnd doesnt change the indices of trials
% plot(timeBall(:,1)- timeStampCa(:,1))


