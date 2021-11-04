function fcoords = arbrot(coords,strike,dip)
%
% ARBROT carries out an arbitrary rotation of points in x,y,z space.
%
%	FCOORDS = ARBROT(COORDS,STRIKE,DIP) uses the x, y, z coordinate
% 	triplets contained in the n x 3 array COORDS and the STRIKE and
%	DIP of the desired rotation (i.e. the vertical axis and horizontal
%	axis rotation, respectively) in ordre to determine the new final
%	coordinates, which are output as coordinate triplets in FCOORDS.
%	To return the points to their original strike, specify [] for 
% 	STRIKE.  To rotate a vertical plane counterclockwise (east dip), 
%	specify a positive DIP, and specify a negative DIP for clockwise
% 	dip rotation.
%
%	The program carries out 5 coordinate transformations:
%	1) After determining the southernmost point, the points are 
%	TRANSLATED so that the southernmost point lies at the origin,
%	2) The points are rotated about their AVERAGE strike so that 
%	they lie in an approximately N-S striking line (plane).
% 	3) The points are rotated about the specified DIP.
%	4) The points are rotated about the specified STRIKE.
%	5) The points are re-translated such that the x,y coordinates
%	of the southernmost point are the same as before the 
%	transformations.
%

warning('off','MATLAB:singularMatrix')

% place the coords into an appropriately shaped array for transformations
coords = [coords'; ones(1, size(coords, 1))];

% find the center point
or = mean(coords, 2);
% create the translation matrix
shift = [[1 0 0;0 1 0;0 0 1;0 0 0] [-or(1:2);0; 1]];

% determine the mean strike of the points
b = [coords(1, :)' ones(size(coords, 2), 1)]\coords(2, :)';
str = 90 - atand(b(1));
north = [cosd(str) -sind(str) 0 0;sind(str) cosd(str) 0 0;0 0 1 0;0 0 0 1];

% create the specified strike matrix
if length(strike) == 0;
	strike = str;
end
strikem = [cosd(strike) sind(strike) 0 0;-sind(strike) cosd(strike) 0 0;0 0 1 0;0 0 0 1];
% dip matrix
di = sign(dip)*(90-abs(dip));
dipm = [cosd(di) 0 -sind(di) 0;0 1 0 0;sind(di) 0 cosd(di) 0;0 0 0 1];
% re-translation matrix
rshift = [[1 0 0;0 1 0;0 0 1;0 0 0] [or(1:2); 0; 1]];

% make the composite transformation
fcoords = rshift*strikem*dipm*north*shift*coords;
fcoords = fcoords(1:3,:)';

warning('on','MATLAB:singularMatrix')