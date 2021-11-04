function [b, uz] = utmconv(a, fi, uz)
% UTMCONV  General UTM conversion function.
%   B = UTMCONV(A, TYPE, UZ) converts the n-by-2 array A from geographic
%   to Cartesian coordinates, or vice versa, using the UTM projection.
%   TYPE specifies the conversion from (long., lat.) to UTM (type = 0) or
%   from UTM to (long., lat.) (type = 1).  UZ specifies the UTM zone, which 
%   should be set as the returned string from UTMZONE for a TYPE = 1
%   conversion and can be either specified or omitted for a TYPE = 0 conversion.
%   If omitted, the zone is determined from the input coordinates.
%

mstruct = defaultm('utm');
if fi == 0 % Lon., lat. to x, y
   if ~exist('uz', 'var');
      mstruct.zone = utmzone(a(:, 2), a(:, 1));
   else
      mstruct.zone = uz;
   end
   mstruct = defaultm(mstruct);
   [x, y] = mfwdtran(mstruct, a(:, 2), a(:, 1));
   b = [x(:) y(:)];
   uz = mstruct.zone;
else
   mstruct.zone = uz;
   mstruct = defaultm(mstruct);
   [lat, lon] = minvtran(mstruct, a(:, 1), a(:, 2));
   b = [lon(:) lat(:)];
end