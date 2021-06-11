function sigSegSuite2p(folder2process,options)
%Runs s2p and/or saves signal, segment and cellsStats files
% folder2process = folder where the .sbx are
% forceS2p : optionnal name-value argument, yes or no
% isCellthresh : optionnal name-value argument, thresholding the iscell stat for keeping cells


%argument parsing new style
arguments
    folder2process (1,:) {mustBeText} 
    options.forceS2p (1,:) string  = "No"
    options.isCellThresh (1,1) {mustBeNumeric} = 0.3
end


%check if the provided path ends with a \ and remove it if yes
if strcmp(folder2process(end),'\')
    folder2process = folder2process(1:end-1);
end

%check if suite2p has been run before
suite2pResultsPath = dir([folder2process '\**\Fall.mat']);
if isempty(suite2pResultsPath) || strcmp(options.forceS2p,'Yes') % if there is no Fall file or if we force S2P
    sbxsuite2p('D:\B10\111020\');
end

%get the Fall file
suite2pResultsPath = dir([folder2process '\**\Fall.mat']);
suite2pResultsPath = [suite2pResultsPath.folder '\' suite2pResultsPath.name ];
s2pRes = load(suite2pResultsPath);

%create the .signals and .segment files
mask  = zeros(s2pRes.ops.Ly,s2pRes.ops.Lx);%create a blank image
goodCells = find(s2pRes.iscell(:,2)>options.isCellThresh); 


for n = 1:length(goodCells)%for each cell
        
        ypix        = s2pRes.stat{goodCells(n)}.ypix(~s2pRes.stat{goodCells(n)}.overlap);%get y positions
        xpix        = s2pRes.stat{goodCells(n)}.xpix(~s2pRes.stat{goodCells(n)}.overlap);%get x positions
        cellSize(n) = length(ypix);
        vert{n} = [xpix' ypix' ]; % /!\ warning :Scanbox segment files has the vertices only, suite2p gives the pixels coordinates. The signal puller uses only the mask anyway
        
        for p = 1:length(ypix) %for each pixel (wtf matlab)
            mask((ypix(p)),xpix(p)) = n;
        end
    
end

%cutting the F and Fneu into separate movies
moviesBoundaries = [0 s2pRes.ops.nframes_per_folder];
for nMovie = 1:length(s2pRes.ops.nframes_per_folder)  
    sigs{nMovie}=s2pRes.F(goodCells,sum(moviesBoundaries(1:nMovie))+1 : ...
        sum(moviesBoundaries(1:nMovie)) + moviesBoundaries(nMovie+1) );
    sigs{nMovie} = sigs{nMovie}';
    sigsBG{nMovie}=s2pRes.Fneu(goodCells,sum(moviesBoundaries(1:nMovie))+1 : ...
        sum(moviesBoundaries(1:nMovie)) + moviesBoundaries(nMovie+1) );
    sigsBG{nMovie}=sigsBG{nMovie}';
end

%saving
files2saveTmp = dir([folder2process '\**\*.sbx']);
files2saveTmp = {files2saveTmp.name};
for f = 1:length(files2saveTmp)
    
    file2save    = files2saveTmp{f};
    file2saveSeg = strrep(file2save,'.sbx','_S2P.segment');
    save([folder2process '\' file2saveSeg ],'mask','vert');
    
    sig          = sigs{f};
    file2saveSig = strrep(file2save,'.sbx','_S2P_signals');
    save([folder2process '\' file2saveSig ],'sig');
    
    sig            = sigsBG{f};
    file2saveSigBG = strrep(file2save,'.sbx','BG_S2P_signals');
    save([folder2process '\' file2saveSigBG ],'sig');
        
end

cellsStats.iscell = s2pRes.iscell(goodCells,2);
cellsStats.cellSize = cellSize; 
file2saveCellsStats = strrep(file2save,'.sbx','CellsStats');
save([folder2process '\' file2saveCellsStats ],'cellsStats');
