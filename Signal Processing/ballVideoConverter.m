function ballVideoConverter(expDir)
cd expDir
files2read = dir('*ball.mj2');

for f = 1:numel(files2read)
    input = [expDir files2read(f).name];
    output = strrep(input,'mj2','mp4');
    command = sprintf('ffmpeg -i %s -vf "scale=trunc(iw/6)*2:trunc(ih/6)*2" %s',input,output);
system(command)
disp('mj2 video converted to mp4')
vidObj = VideoReader(output);
f=0;
while hasFrame(vidObj)
    f=f+1;
    vidFrame = readFrame(vidObj);
        vidFrame = vidFrame(:,:,1) ;

    data(:,:,f) = vidFrame;
end
time = 0:1/vidObj.FrameRate:vidObj.Duration-1/vidObj.FrameRate;  

save(strrep(output,'mp4','mat'),'data',"time")
clear vidObj
end
  
