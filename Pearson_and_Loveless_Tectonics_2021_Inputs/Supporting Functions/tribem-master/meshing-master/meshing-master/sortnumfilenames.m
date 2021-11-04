function sfiles = sortnumfilenames(files)
% SORTNUMFILENAMES  Sorts filenames in true numerical order.
%   SORTFNUMFILENAMES(FILES) sorts a structure of filenames returned from
%   DIR, given as FILES, in true numerical order (i.e., file1.txt, file2.txt,...) 
%   rather than the default result from DIR of file1.txt, file10.txt,...
%
%   Any filenames containing a number will be properly sorted. Filenames
%   containing no digits will be placed at the end of the list but will 
%   remain in the same order relative to each other as they appear in FILES.
%
%   SFILES = SORTNUMFILENAMES(FILES) returns the sorted file structure to SFILES.
%

nf = numel(files); % Number of files
N = zeros(nf, 1); % Allocate space for sorting matrix
for i = 1:nf
   ns = regexp(files(i).name, '\d+', 'match'); % Check for number as string
   if ~isempty(ns) % If a number exists, 
      N(i) = str2num(ns{:}); % Convert it to numeric
   else % If it's entirely text, 
      N(i) = 1e10+i; % Give this file a large value so it appears at the end of the list
   end
end
[~, Ns] = sort(N); % Find indices of sorted numerical values. Text-only names are placed last, but in the original order relative to each other
sfiles = files(Ns); % Sort files
