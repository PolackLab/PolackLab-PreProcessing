function sbxsuite2p(fname)
 
str = [];
 
if iscell(fname)
    for i=1:length(fname)
        str = [str ' ' getpath(fname{i})];
    end
else
    str = getpath(fname);
end
 
str = strrep(str,'\','/');
pyscript = strrep(which('suite2p_process.py'),'\','/');     % find the processing script
cmd = sprintf('CALL conda.bat activate suite2p && python.exe %s %s & ', pyscript, str);
system(cmd);
 
 
function p = getpath(fn)
 
try
    s = what(fn);
    p = s.path;
catch
    error('No such file');
end