function varargout = contour2xyz(c)
% contour2xyz  Converts contour array to x, y, z
%   [x, y, z] = contour2xyz(c) converts the contour array c
%   to vectors x, y, and z.
% 
%   xyz = contour2xyz(c) converts the contour array to the  
%   3-column array xyz.
%   
%

% Allocate space for arrays
[x, y, z] = deal(zeros(size(c, 2), 1));

i = 0; % Counting variable along columns of contours
j = 1; % Count number of contour lines
while i < size(c, 2)
   cn = c(2, i+1); % Get length of contour
   dep = c(1, i+1); % Get contour interval
   x(i+(1:cn)) = c(1, i+(2:cn+1)); % Insert x coordinates
   y(i+(1:cn)) = c(2, i+(2:cn+1)); % Insert y coordinates
   z(i+(1:cn)) = dep; % Replicate z coordinate
   i = i+1+cn; % Increment counter
   j = j+1; % Increment number of contour lines
end
keep = x ~= 0;
x = x(keep);
y = y(keep);
z = z(keep);

if nargout == 1
   varargout{1} = [x(:), y(:), z(:)];
elseif nargout == 3
   varargout{1} = x;
   varargout{2} = y;
   varargout{3} = z;   
end
