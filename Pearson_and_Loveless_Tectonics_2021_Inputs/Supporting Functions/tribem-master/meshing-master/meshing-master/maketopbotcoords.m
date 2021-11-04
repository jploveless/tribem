function maketopbotcoords(file, depth)
% maketopbotcoords   Makes files for top and bottom coordinates from an ArcGIS text file.
%   maketopbotcoords(FILE, DEPTH) makes folders containing x, y, z coordinates for 
%   the top of a fault and the bottom (with uniform z coordinate specified by input
%   arugment DEPTH), based on a FILE containing surface trace coordinates exported from
%   ArcGIS. Individual text files will be generated for distinct features within the 
%   coordinate file. Folders will be generated to hold the top and bottom coordiante files,
%   named FILE_topcoords and FILE_botcoords, respectively. 
%
%   Use the function rotatefaultcoords to assign dips to each fault (by changing the
%   bottom coordinates) and function gmshfaults to create triangulated meshes of all
%   faults.
%

% ArcGIS file style 1: Feature coordinates start with FID and end with "END"
a = opentxt(file);
ends = strfind(a(:, 1)', 'E') - 1;
if ~isempty(ends)
   astyle = 1;
   if diff(ends(end), ends(end-1)) == 1
      ends = ends(1:end-1);
   end
   begs = [2 ends(1:end-1) + 3];
else
% ArcGIS file style 2: Each line contains feature x, y, FID
   astyle = 2;
   if a(1) == 'X'
      fl = 2;
   else
      fl = 1;
   end
   a = str2num(a(fl:end, :));
   d3 = diff(a(:, 3));
   ends = [find(d3==1); size(a, 1)];
   begs = [1; ends(1:end-1)+1];
end

topdir = [file(1:end-4), '_topcoords'];
botdir = [file(1:end-4), '_botcoords'];
if ~exist(topdir, 'dir')
   mkdir(topdir)
end

if ~exist(botdir, 'dir')
   mkdir(botdir)
end

for i = 1:numel(begs)
   if astyle == 1
      tout = str2num(a(begs(i):ends(i), :));
   else
      tout = [a(begs(i):ends(i), 1:2), zeros(ends(i)-begs(i)+1, 1)];
   end

   tout = unique(tout, 'rows', 'stable');
   % Try to sort coordinates along a single direction
%   [~, di] = sort(mag(tout, 2));
%   tout = tout(di, :);

   bout = [tout(:, 1:2) -abs(depth)*ones(size(tout, 1), 1)];
   save(sprintf('%s/f%g.txt', topdir, i), 'tout', '-ascii');
   save(sprintf('%s/f%g.txt', botdir, i), 'bout', '-ascii');
end