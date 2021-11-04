function p = gmshsurface(edgex, edgey, edgez, top, bot, s1, s2, sz)
% gmshplane  Makes a triangulated nonplanar fault in Gmsh. 
%    p = gmshsurface(edgex, edgey, edgez, top, bot, s1, s2, sz) creates a 
%    triangulated mesh structure using Gmsh for the surface defined
%    by the edge coordiantes input variables. edgex, edgey, edgez are 
%    vectors giving all edge coordinates. top, bot, s1, and s2 are 
%    vectors giving the indices of the edge coordinate arrays to define
%    the top, bottom, and sides of the surface. The element size is given 
%    by sz. 
%   

np = length(edgex); % Number of points
alledge = [top(:); s1(:); bot(:); s2(:)]; % All edges
edgeend = cumsum([length(top), length(s1), length(bot), length(s2)]); % Edge indices
edgebeg = [1 edgeend(1:end-1)]; 

% Write coordinates to Gmsh .geo file
fn = sprintf('junk%g', rem(now, 1));
fid = fopen([fn '.geo'], 'w');
fprintf(fid, 'charl = %g;\n', sz); % Write element size
% Write points: top then bottom
fprintf(fid, 'Point(%g) = {%g, %g, %g, charl};\n', [1:np; edgex(:)'; edgey(:)'; edgez(:)']);
% Write lines: top, side 1, bottom, side 2
for i = 1:4,
   edge = alledge(edgebeg(i):edgeend(i));
   fprintf(fid, 'BSpline(%g) = {%g', i, edge(1));
   fprintf(fid, ', %g', edge(2:end));
   fprintf(fid, '};\n');
end
% Circulate lines and create plane
fprintf(fid, 'Line Loop(1) = {1, 2, 3, 4};\nSurface(1) = {1};\n');
fclose(fid);

% Mesh using Gmsh
p = gmsh([fn '.geo']);

% Remove temp files
system(sprintf('rm %s*', fn));
