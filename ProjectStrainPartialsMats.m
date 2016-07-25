function Gp = ProjectStrainPartials(G, strike, dip, rake);
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
dip(dip > 90) = dip(dip > 90) - 90;

% Check rake
if ~exist('rake', 'var')
   rake = zeros(size(dip));
end

% Trig. function arrays
sins = repmat(sind(strike), 1, size(G, 2));
coss = repmat(cosd(strike), 1, size(G, 2));
sined = repmat(sind(dip), 1, size(G, 2));
cosed = repmat(cosd(dip), 1, size(G, 2));
cosr = repmat(cosd(rake), 1, size(G, 2));
sinr = repmat(sind(rake), 1, size(G, 2));

% Extract tensor components
uxx = G(1:6:end, :);
uyy = G(2:6:end, :);
uzz = G(3:6:end, :);
uxy = G(4:6:end, :);
uxz = G(5:6:end, :);
uyz = G(6:6:end, :);

% Do the rotation, R'*D'*S'*T*S*D*R (used symbolic toolbox)
ud = -sined.*(coss.*(uxz.*(sined) - uxx.*(coss).*(cosed) + uxy.*(cosed).*(sins)) - sins.*(uyz.*(sined) - uxy.*(coss).*(cosed) + uyy.*(cosed).*(sins))) - cosed.*(uzz.*(sined) - uxz.*(coss).*(cosed) + uyz.*(cosed).*(sins));
us = sined.*(coss.*(uxy.*(coss) + uxx.*(sins)) - sins.*(uyy.*(coss) + uxy.*(sins))) + cosed.*(uyz.*(coss) + uxz.*(sins));
un = sined.*(coss.*(uxz.*(cosed) + uxx.*(coss).*(sined) - uxy.*(sins).*(sined)) - sins.*(uyz.*(cosed) + uxy.*(coss).*(sined) - uyy.*(sins).*(sined))) + cosed.*(uzz.*(cosed) + uxz.*(coss).*(sined) - uyz.*(sins).*(sined));
 
 
% Write the projected partials
Gp = zeros(size(G, 1)/2, size(G, 2));
Gp(1:3:end, :) = us;
Gp(2:3:end, :) = ud;
Gp(3:3:end, :) = un;

