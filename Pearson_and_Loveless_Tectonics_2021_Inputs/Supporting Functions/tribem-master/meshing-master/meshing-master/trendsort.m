function sorted = trendsort(coords, thresh)
% trendsort  Sorts x, y coordinates along a similar trend direction
%   sorted = trendsort(coords) sorts the n-by-2 array coords, containing 
%   x, y coordinate pairs so that they follow a single general trend (i.e.,
%   east to west, or south to north).
%

% First sort based on distance from a common origin
[~, di] = sort(mag(coords, 2));
sorted = coords(di, :);

% Calculate trends of sorted coordinate pairs
dc = diff(sorted);
trends = atan2d(dc(:, 1), dc(:, 1));

