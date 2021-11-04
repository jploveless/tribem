function p = gmsh(geofile)
% gmsh  Calls gmsh to mesh a 2-D surface.
%   p = gmsh(geofile) will call Gmsh to mesh a 2-D fault
%   surface from the geometry file geofile and return the
%   fault to structure p. 
%

% Make sure the .geo file has .geo extension
[p, f, e] = fileparts(geofile);
if isempty(e)
   geofile = [geofile '.geo'];
else
   if ~startsWith(e, 'geo')
      geofile = [p, f, '.geo'];
   end
end

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

% Do the meshing
system(sprintf('%s/gmsh -2 %s -o %s.msh -v 0 > junk', gmshpath, geofile, geofile(1:end-4)));

% Read the mesh
p = ReadPatches(sprintf('%s.msh', geofile(1:end-4)));

% Remove temp files
system('rm junk');
