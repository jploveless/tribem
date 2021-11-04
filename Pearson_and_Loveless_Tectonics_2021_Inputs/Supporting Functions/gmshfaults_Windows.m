function allstring = gmshfaults_Windows(topdir, botdir, sz, path)
% gmshfaults  Meshes faults in 3D using Gmsh.
%   gmshfaults(BOTDIR, ELSIZE) meshes all faults whose bottom coordinates
%   are contained in text files within the folder BOTDIR, using the 
%   element size specified using ELSIZE.  Geometry (.geo) and mesh (.msh)
%   files are saved to the bottom coordinate directory.
%
%   ALL = gmshfaults(...) returns a string to the command line that can be
%   used to read all faults in using ReadPatches.m, i.e.
%   >> p = ReadPatches(ALL);
% 
%

% Hard code full path to top coordinates directory here
%topdir = ['~' filesep, 'Desktop', filesep,'chile', filesep,'Manual_Faultloc_All', filesep,'Manual_Faultloc_Originals', filesep,'Manual_Faultloc_16km'];

% List all bottom coordinate files
files = dir([botdir filesep '*.txt']);

for i = 1:size(files, 1)
   % Read bottom coordinates
   botcoords = load([botdir filesep files(i).name]);
   % Get number of coordinates
   nc = size(botcoords, 1);
   % Read top coordinates
   topcoords = load([topdir filesep files(i).name]);
   % Write the Gmsh geometry file
   geofile = [botdir filesep files(i).name(1:end-4) '.geo'];
   fid = fopen(geofile, 'w');
   fprintf(fid, 'charl = %g;\n', sz);
   fprintf(fid, 'Point(%g) = {%g, %g, %g, charl};\n', [1:2*nc; [topcoords(:, 1); flipud(botcoords(:, 1))]'; [topcoords(:, 2); flipud(botcoords(:, 2))]'; [topcoords(:, 3); flipud(botcoords(:, 3))]']);
   fprintf(fid, 'CatmullRom(%g) = {%g:%g};\n', [1:2; [1 nc+1; nc 2*nc]]);
   fprintf(fid, 'CatmullRom(%g) = {%g,%g};\n', [3:4; [nc 2*nc; nc+1 1]]);
   fprintf(fid, 'Line Loop(1) = {1, 3, 2, 4};\nRuled Surface(1) = {1};\n');
   fclose(fid);
   fid = fopen('infile.txt','rt') ;

   % Mesh using Gmsh (default)
   % system(sprintf('/Applications/Gmsh/Gmsh.app/Contents/MacOS/gmsh -2 %s -o %s.msh -v 0 > junk', geofile, geofile(1:end-4)));
   % Mesh using Gmsh (Linux)
%    system(sprintf('gmsh -2 %s -o %s.msh -v 0 > junk', geofile, geofile(1:end-4)));
    % Mesh using Gmsh (Windows)
       system(sprintf(strcat(path,' -2 %s -o %s.msh -v 0 > junk'), geofile, geofile(1:end-4)));

end

% Make a space-separated string of all .msh files
alldir = dir([botdir '/*.msh']); 
xlsfiles={alldir.name};
[~,idx]=natsort(xlsfiles);
alldir=alldir(idx);
allstring = sprintf('%s ', alldir.name); 
allstring = allstring(1:end-1); % Remove trailing space