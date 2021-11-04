function contour2geo(lon, lat, dep, esize, file, gvers)
% contour2geo  Converts contour lines to a .geo file for Gmsh.
%   contour2geo(lon, lat, dep, esize, file) converts the contour
%   lines whose coordinates are defined by vectors lon, lat, dep
%   into a .geo file for meshing using Gmsh. esize specifies the 
%   element size, and file gives the filename with full path. 
%
%   The lon, lat of contour lines are assumed to be ordered, 
%   and contour lines are separated based on changes in depth.
%
%   The geographic coordinates are converted to geocentric 
%   coordinates before writing to the file. These can be 
%   converted back to lon, lat, dep using mshxyz2geo.m. 

% Identify breaks in contour lines based on depth
dd = diff(dep);
maxdd = max(abs(dd));
breaks = find(abs(dd) > 0.2*maxdd);
begs = [1; breaks+1];
ends = [breaks; size(dep, 1)];

% Check for Gmsh version ID
if exist('gvers', 'var')
   if gvers >= 4
      surfname = 'Surface';
   else
      surfname = 'Ruled Surface';
   end
else
   surfname = 'Surface'; % Assume newer version
end

% Convert all coordinates to geocentric Cartesian
[x, y, z] = long_lat_to_xyz(deg2rad(lon), deg2rad(lat), 6371-abs(dep));

% Prep file for writing
fid = fopen([file '.geo'], 'w');
fprintf(fid, 'esize = %g;\n', esize);
% Write .geo points
fprintf(fid, 'Point(%g) = {%g, %g, %g, esize};\n', [1:size(lon, 1); x'; y'; z']);
% Write contour lines
for i = 1:length(begs)
   fprintf(fid, 'CatmullRom(%g) = {%g:%g};\n', i, begs(i), ends(i));
end
j = i;
% Write end lines, circulations, and surfaces
for i = 1:length(begs)-1
   % Depth connection 1
   fprintf(fid, 'CatmullRom(%g) = {%g,%g};\n', j+(2*i-1), begs(i), begs(i+1));
   % Depth connection 2
   fprintf(fid, 'CatmullRom(%g) = {%g,%g};\n', j+(2*i-0), ends(i), ends(i+1));
   % Line loop
   fprintf(fid, 'Line Loop(%g) = {%g, %g, -%g, -%g};\n', i, i, j+(2*i-0), i+1, j+(2*i-1));
   % Surface
   fprintf(fid, '%s(%g) = {%g};\n', surfname, i, i);
end

fclose(fid);