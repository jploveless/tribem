function mnames = mshfilenames(dir)
% mshfilenames  Makes a space-separated string of a directory's .msh filenames
%   mnames = mshfilenames(dir) creates a space-delimited string of the .msh
%   filenames within folder dir. mnames is a string suitable for use with 
%   ReadPatches.m.
%

mnames = '';
files = dir('*.msh')'

for i = 1:numel(files)
   mnames = sprintf('%s %s/%s', name, files(i).folder, files(i).name);
end

mnames = mnames(2:end);