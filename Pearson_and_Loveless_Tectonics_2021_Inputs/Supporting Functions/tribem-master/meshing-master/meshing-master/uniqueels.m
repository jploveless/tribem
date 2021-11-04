function pu = uniqueels(p)
% uniqueels   Cleans a patch structure retaining unique elements.
%   pu = uniqueels(p) cleans patch structure p leaving only unique
%   elements, returning them to structure pu. 

% Get node coordinates
px = unstack3(p.c(p.v, 1));
py = unstack3(p.c(p.v, 2));
pz = unstack3(p.c(p.v, 3));

% Check unique nodes
[uc, ic] = unique([px, py, pz], 'rows');

% Extract unique structure
pu = patchsubset(p, ic);

