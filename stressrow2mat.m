function mat = stressrow2mat(row)
% stressrow2mat  Converts stress tensor components from rows to matrices.
%   mat = stressrow2mat(row) converts stress tensor components from an
%   n-by-6 matrix where components are arranged in rows to a 3n-by-3 matrix
%   where components are arranged in 3-by-3 blocks. 
%
%   The tensor components in the input array are assumed to be ordered
%   [xx, yy, zz, xy, xz, yz]
%

% Number of tensors
nrows = size(row, 1);

% Allocate space for reshaped matrix
mat = zeros(3*nrows, 3);

% Insert row entries in correct row, column locations
mat(1:3:end, :) = row(:, [1 4 5]);
mat(2:3:end, :) = row(:, [4 2 6]);
mat(3:3:end, :) = row(:, [5 6 3]);
