function p = contour2p(lon, lat, dep, esize)
% contour2p  Converts contour lines to mesh structure.
%   p = contour2p(lon, lat, depth, esize, geofile) converts the
%   ordered contour lines with coordinates lon, lat, depth into a 
%   mesh structure p, with elements of size esize, meshed using
%   Gmsh. 
% 

% Check for preferences file
if exist('gmshfaultspref.mat', 'file') ~= 0 % If this .mat file exists, 
   load('gmshfaultspref.mat', 'gmshpath') % Load it
else % If not, 
   if ismac
      if exist('/Applications/Gmsh.app/Contents/MacOS/gmsh', 'file') % Check for default install location
         gmshpath = '/Applications/Gmsh.app/Contents/MacOS/';
      else
         gmshpath = ''; % Or ask for install location
         while ~exist([gmshpath filesep 'gmsh'], 'file')
            gmshpath = input('Enter path to Gmsh application: ');
         end
      end
      % Save Gmsh path to preferences file, to be read in future runs
      gmfp = fileparts(which('gmshfaults'));
      save([gmfp filesep 'gmshfaultspref.mat'], 'gmshpath');
   elseif ispc || (isunix && ~ismac)
      gmshpath = ''; % Or ask for install location
      while ~exist([gmshpath filesep 'gmsh.exe'], 'file')
         gmshpath = input('Enter path to Gmsh application: ');
      end
      % Save Gmsh path to preferences file, to be read in future runs
      gmfp = fileparts(which('gmshfaults'));
      save([gmfp filesep 'gmshfaultspref.mat'], 'gmshpath');
   end
end

% Write geo file, first checking version
geofile = 'from_contour2p';
[~, gvers] = system(sprintf('%s/gmsh -version', gmshpath));
gvers = regexp(gvers, '\d+.\d+.\d+', 'match');
gvers = str2num(gvers{end}(1));
contour2geo(lon, lat, dep, esize, geofile, gvers)

% Do the meshing
system(sprintf('%s/gmsh -2 %s.geo -o %s.msh -v 0 > junk', gmshpath, geofile, geofile));

% Read in the mesh coordinates
p = ReadPatches([geofile '.msh']);

% Convert from geocentric to geographic
[lon, lat] = xyz_to_long_lat(p.c(:, 1), p.c(:, 2), p.c(:, 3)); % Long, lat
cdis = mag(p.c, 2); % Distance from earth center
p.c = [rad2deg(lon), rad2deg(lat), cdis-6371];