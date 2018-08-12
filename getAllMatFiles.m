function [filelist, nfile] = getAllMatFiles(foldername, ext)
% get the number and name of all files with a specific extension in a folder
d = dir([foldername '\*.' ext ]);
nfile = size(d,1);
filelist = {d.name};
end
