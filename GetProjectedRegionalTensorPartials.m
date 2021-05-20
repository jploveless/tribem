function Gp = ProjectStrainPartials(strike, dip, rake);
% PROJECTSTRAINPARTIALS  Projects strain partial components onto a fault.
%
%   PROJECTSTRAINPARTIALS(G, S, D) projects the matrix of strain partials G onto
%   the fault surfaces defined by strike, S, and dip, D, arrays.  G is a 6.*n-by-3.*n
%   array relating the 6 strain components to the 3 slip components of n faults, and
%   S and D are n-by-1 arrays defining each fault's geometry.
%
%   GP = PROJECTSTRAINPARTIALS(G, S, D) returns the projected matrix to GP, which is a 
%   3.*n-by-3.*n array with rows describing traction in the strike direction, in the dip
%   direction, and in the fault-normal direction.
%
%   GP = PROJECTSTRAINPARTIALS(G, S, D, R) also uses a rake vector R and, instead of 
%   returning traction in the strike and dip direction, returns traction in the rake
%   direction and the rake perpendicular direction.
%

% Handle dips > 90
strike(dip > 90) = wrapTo360(strike(dip > 90) + 180);
dip(dip > 90) = 180 - dip(dip > 90);

% Check rake
if ~exist('rake', 'var')
   rake = zeros(size(dip));
end

% Trig. function arrays
sins  = sind(strike);
coss  = cosd(strike);
sined = sind(dip);
cosed = cosd(dip);
cosr  = cosd(rake);
sinr  = sind(rake);

% Allocate space for partials matrix
Gp = zeros(3*numel(strike), 6);

% Do the rotation, R'*D'*S'*T*S*D*R (used symbolic toolbox)
Gp(1:3:end, :) = [sined.*coss.*sins, -sined.*sins.*coss, 0*sins, sined.*(coss.*coss - sins.*sins), cosed.*sins, cosed.*coss];
Gp(2:3:end, :) = -[sined.*coss.*coss.*cosed, sined.*sins.*cosed.*sins, -cosed.*sined, -sined.*(coss.*cosed.*sins + sins.*coss.*cosed), -sined.*coss.*sined + cosed.*coss.*cosed, sined.*sins.*sined - cosed.*cosed.*sins];
Gp(3:3:end, :) = [sined.*coss.*coss.*sined, sined.*sins.*sins.*sined, cosed.*cosed, sined.*(-coss.*sins.*sined - sins.*coss.*sined), sined.*coss.*cosed + cosed.*coss.*sined, -sined.*sins.*cosed - cosed.*sins.*sined];
