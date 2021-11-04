function p = gmshlistric(edgex, edgey, edgez, top, bot, ls1, ls2, fs1, fs2, sz)
% gmshplane  Makes a triangulated listric fault in Gmsh. 
%    p = gmshsurface(edgex, edgey, edgez, top, bot, ls1, ls2, fs1, fs2, sz)
%    creates a triangulated mesh structure using Gmsh for the surface defined
%    by the edge coordiantes input variables. edgex, edgey, edgez are 
%    vectors giving all edge coordinates. top, bot, ls1, ls2, fs1, and fs2 are 
%    vectors giving the indices of the edge coordinate arrays to define
%    the top, bottom, listric sides, and flat sides of the surface. The element
%    size is given by sz. 
%   

np = length(edgex); % Number of points
nrows = size(ls1, 2)/2; % Number of listric rows
alledge = [top(:); ls1(:); fs1(:); bot(:); fs2(:); ls2(:)]; % All edges
nedge = size(alledge, 1)/2; % Number of edges

% Write coordinates to Gmsh .geo file
%fn = sprintf('junk%g', rem(now, 1));
fn = 'junk';
fid = fopen([fn '.geo'], 'w');
charl = sz*ones(1, np); % Create vector of element size
charl(bot) = 10*sz; % Make nodes along bottom have a larger size
% Write points: top then bottom
fprintf(fid, 'Point(%g) = {%g, %g, %g, %g};\n', [1:np; edgex(:)'; edgey(:)'; edgez(:)'; charl]);
% Write lines: top, listric side 1, flat side 1, bottom, flat side 2, listric side 2
for i = 1:nedge
   fprintf(fid, 'Line(%g) = {%g, %g};\n', i, alledge(2*i-1), alledge(2*i));
end

% Extra lines that separate panels
for i = 1:nrows
   fprintf(fid, 'Line(%g) = {%g, %g};\n', nedge+i, ls1(length(ls1)-2*(i-1)), ls2(2*i-1));
end

% Circulate lines and create surfaces
fprintf(fid, 'Line Loop(1) = {%g, %g, %g, %g};\nSurface(1) = {1};\n', [nrows+(2:4), -(2*nrows+5)]); % Flat part of fault
for i = 1:nrows % Listric panels
   if i == 1
      fprintf(fid, 'Line Loop(%g) = {%g, %g, %g, %g};\nSurface(%g) = {%g};\n', 1+i, 2, nedge+nrows, nedge, 1, 1+i, 1+i);
   else 
      fprintf(fid, 'Line Loop(%g) = {%g, %g, %g, %g};\nSurface(%g) = {%g};\n', 1+i, nrows+3-i, 2*nrows+i+3, nrows+i+3, -(2*nrows+i+4), 1+i, 1+i);
   end
end
fclose(fid);

% Mesh using Gmsh
p = gmsh([fn '.geo']);

% Remove temp files
%system(sprintf('rm %s*', fn));
