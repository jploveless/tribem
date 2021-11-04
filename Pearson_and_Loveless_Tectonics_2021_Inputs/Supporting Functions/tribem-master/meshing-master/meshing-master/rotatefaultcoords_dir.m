function rotcoords=rotatefaultcoords(dirname, dip)
% rotatefaultcoords  Rotates a set of fault coordinates to a specified dip.
%
%   rotatefaultcoords(DIRNAME, DIP) rotates the coordinates contained in
%   all text files within directory DIRNAME to represent a fault of a 
%   specified DIP.  Negative DIP values give westward dips, positive give
%   eastward dips.  The text files are assumed to contain X, Y, and Z 
%   coordinates of the trace of a vertical fault at depth.  The rotated 
%   coordinates will be saved to a text file within a newly created 
%   directory with "_DIPdip" appended.
%

% Create a new directory into which rotated coordinate files will be placed
rotdirname = sprintf('%s_truedip', dirname);
if ~exist(rotdirname, 'dir')
   mkdir(rotdirname);
end
% Get a listing of all files in the specified directory
files = dir([dirname filesep '*.txt']);
% Sort files numerically
files = sortnumfilenames(files);

% Expand scalar dip to one value per file
if numel(dip) == 1
   dip = dip.*ones(numel(files), 1);
end

% For each file...
for i = 1:size(files, 1)
   % Load the coordinates
   coords = load([dirname filesep files(i).name]);
   % Rotate the coordinates by the specified dip
   rotcoords = arbrot(coords, [], dip(i));
   % Save the rotated coordinates to the newly created directory
   save([rotdirname filesep files(i).name], 'rotcoords', '-ascii');
end
   