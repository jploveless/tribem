function p = gmshplane(tx, ty, tz, maxz, dip, sz)
% gmshplane  Makes a triangulated planar fault in Gmsh. 
%    p = gmshplane(tx, ty, tz, L, W, strike, dip, size) creates a 
%    triangulated mesh structure using Gmsh for the plane defined
%    by the input variables. tx, ty, tz are 2-element vectors giving
%    the coordinates of the end points of the tip line; maxz gives 
%    the maximum depth; dip gives the dip, with positive to the right 
%    looking from endpoint 1 to endpoint 2; and size gives the nominal
%    element size ("characteristic length") in Gmsh. 
%   

% Determine deep coordinates based on dip and maximum depth
strike = atan2d(diff(tx), diff(ty));
horzd = (abs(maxz) - abs(tz))./tand(dip);
dx = tx + cosd(strike).*horzd;
dy = ty + sind(strike).*horzd;
dz = [maxz; maxz];

% Write coordinates to Gmsh .geo file
fid = fopen('junk.geo', 'w');
fprintf(fid, 'charl%g = %g;\n', [1:length(sz); sz(:)']); % Write element size
szi = [1 1 length(sz) length(sz)]; % Element size index array
% Write points: top then bottom
fprintf(fid, 'Point(%g) = {%g, %g, %g, charl%g};\n', [1:4; [tx(:)', dx(:)']; [ty(:)', dy(:)']; [tz(:)', dz(:)']; szi]);
% Write lines: top, side 1, bottom, side 2
fprintf(fid, 'Line(%g) = {%g,%g};\n', [1, 1, 2, 2, 2, 4, 3, 4, 3, 4, 3, 1]);
% Circulate lines and create plane
fprintf(fid, 'Line Loop(1) = {1, 2, 3, 4};\nRuled Surface(1) = {1};\n');
fclose(fid);

% Mesh using Gmsh
p = gmsh('junk.geo');

% Remove temp files
system('rm junk*');
