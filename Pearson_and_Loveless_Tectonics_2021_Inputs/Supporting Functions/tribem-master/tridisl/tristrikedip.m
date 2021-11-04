function [s, d, nv] = tristrikedip(p,t)

% TRISTRIKEDIP: Strike and dip of triangles assuming counter-clockwise
% (CCW) node ordering.
%
%  P  : Nx3 array of XYZ node co-ordinates
%  T  : Mx3 array of triangles as indices into P
%  S  : Mx1 array of triangle strikes
%  D  : Mx1 array of triangle dips

d12 = p(t(:,2),:)-p(t(:,1),:);
d13 = p(t(:,3),:)-p(t(:,1),:);
nv = cross(d12, d13, 2);
nv = nv./repmat(mag(nv, 2), 1, 3);
[s, d] = cart2sph(nv(:, 1), nv(:, 2), nv(:, 3));
s = wrapTo360(-rad2deg(s));
d = 90 - rad2deg(d);

end      % tristrikedip()